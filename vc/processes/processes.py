from vc import *

# 仅计算细胞每次接触时产生的trimers，并不产生PD效果
class ComputeTrimersProcess:
  # contact_duration：T细胞和靶细胞接触时间
  def __init__(self, system, effector, affected, drug, compartments,
               contact_duration = 10*units.MIN, contact_area = 1*units.um**2,
               time_threshold = 30*units.s):
    system.reports[self] = []
    
    self.system = system
    self.effector = effector
    self.affected = affected
    self.drug = drug
    self.compartments = compartments
    
    self.contact_duration = contact_duration.number(units.h)
    self.contact_area = contact_area.number(units.um**2)
    self.time_threshold = time_threshold.number(units.h)
    
    _effector = system.cells.index(effector)
    _affected = system.cells.index(affected)
    _compartments = [system.compartments.index(compartment) for compartment in compartments]
    names_effector, names_affected, _mat_on, _mat_off = drug._get_trans_mats(effector, affected, system.solutes)
    analytes_effector = [f"({effector.name}){name_effector}" for name_effector in names_effector]
    analytes_affected = [f"({affected.name}){name_affected}" for name_affected in names_affected]
    _analytes_effector = [system.analyteses_bound[_effector].index(analyte_effector) for analyte_effector in analytes_effector]
    _analytes_affected = [system.analyteses_bound[_affected].index(analyte_affected) for analyte_affected in analytes_affected]
    self.precomp = _effector, _affected, _compartments, _analytes_effector, _analytes_affected, _mat_on, _mat_off
  
  def __call__(self, system, t):
    assert system == self.system, "The system to apply process on is different than the one specified at process instantiation!"
    
    # 此阶段计算p -- 一次接触能形成免疫突触的概率
    effector, affected, compartments, drug = self.effector, self.affected, self.compartments, self.drug
    _effector, _affected, _compartments, _analytes_effector, _analytes_affected, _mat_on, _mat_off = self.precomp
    yy = np.einsum("ik,jk->ijk", system.y[_effector][np.ix_(_analytes_effector, _compartments)], system.y[_affected][np.ix_(_analytes_affected, _compartments)]) # 二维浓度乘积张量，维度是len(analytes_effector) x len(analytes_affected) x len(compartments)
    trimers_formed = np.einsum("ijk,ij->ijk", yy, _mat_on) * self.contact_area * self.contact_duration # 总共产生多少各类型的trimers
    probs_above_threshold = np.exp(-_mat_off * self.time_threshold) # 各类型的trimers能维持超过阈值时间的概率
    trimers_above_threshold = np.einsum("ijk,ij->k", trimers_formed, probs_above_threshold)
    trimers_below_threshold = np.einsum("ijk,ij->k", trimers_formed, 1-probs_above_threshold)
    total_drug_on_effector = system.get_total_drug_on_cell(compartments, effector, drug)
    total_drug_on_affected = system.get_total_drug_on_cell(compartments, affected, drug)
    
    buffer = {"t":system.t, "trimers_above_threshold":trimers_above_threshold, "trimers_below_threshold":trimers_below_threshold, "total_drug_on_effector":total_drug_on_effector, "total_drug_on_affected":total_drug_on_affected, "yy":yy}
    system.reports[self].append(buffer)


class CytotoxicityProcess:
  # contact_freq
  # contact_duration：T细胞和靶细胞接触时间
  # damage_per_hit：每次攻击造成的损伤
  # cv_damage_per_hit：每次攻击造成的损伤的coefficient of variation
  # repair：修复损伤的速度
  def __init__(self, system, effector, affected, drug, compartments,
               n = 100000, #计算hp时的颗粒度
               contact_freq = 1 / (3*units.h) / (5e4/units.ml), # 假设ET比1:1的条件下每3小时接触一次
               contact_duration = 10*units.MIN, contact_area = 1*units.um**2,
               time_threshold = 30*units.s, prob_synapse = None, # 免疫突触质量体现在这个参数
               damage_per_hit = None, cv_damage_per_hit = None, repair = None): #最后几个参数和细胞有关，默认读取细胞的属性，但可以在这里覆盖
    self.system = system
    self.effector = effector
    self.affected = affected
    self.drug = drug
    self.compartments = compartments
    
    self.n = n
    self.contact_freq = contact_freq.number(units.ml/units.h)
    self.contact_duration = contact_duration.number(units.h)
    self.contact_area = contact_area.number(units.um**2)
    self.time_threshold = time_threshold.number(units.h)
    
    self.prob_synapse = prob_synapse if damage_per_hit is not None else drug.synapse_probability
    
    self.damage_per_hit = damage_per_hit if damage_per_hit is not None else affected.damage_per_hit
    self.cv_damage_per_hit = cv_damage_per_hit if cv_damage_per_hit is not None else affected.cv_damage_per_hit
    self.repair = repair.number(1/units.h) if repair is not None else affected.repair.number(1/units.h)
    
    _effector = system.cells.index(effector)
    _affected = system.cells.index(affected)
    _compartments = [system.compartments.index(compartment) for compartment in compartments]
    names_effector, names_affected, _mat_on, _mat_off = drug._get_trans_mats(effector, affected, system.solutes)
    analytes_effector = [f"({effector.name}){name_effector}" for name_effector in names_effector]
    analytes_affected = [f"({affected.name}){name_affected}" for name_affected in names_affected]
    _analytes_effector = [system.analyteses_bound[_effector].index(analyte_effector) for analyte_effector in analytes_effector]
    _analytes_affected = [system.analyteses_bound[_affected].index(analyte_affected) for analyte_affected in analytes_affected]
    self.precomp = _effector, _affected, _compartments, _analytes_effector, _analytes_affected, _mat_on, _mat_off
    self.hp = np.ones((self.n, len(compartments)))
  
  def __call__(self, system, t):
    assert system == self.system, "The system to apply process on is different than the one specified at process instantiation!"
    
    # 此阶段计算p -- 一次接触能形成免疫突触的概率
    _effector, _affected, _compartments, _analytes_effector, _analytes_affected, _mat_on, _mat_off = self.precomp
    yy = np.einsum("ik,jk->ijk", system.y[_effector][np.ix_(_analytes_effector, _compartments)], system.y[_affected][np.ix_(_analytes_affected, _compartments)]) # 二维浓度乘积张量，维度是len(analytes_effector) x len(analytes_affected) x len(compartments)
    trimers_formed = np.einsum("ijk,ij->ijk", yy, _mat_on) * self.contact_area * self.contact_duration # 总共产生多少各类型的trimers
    probs_above_threshold = np.exp(-_mat_off * self.time_threshold) # 各类型的trimers能维持超过阈值时间的概率
    trimers_above_threshold = np.einsum("ijk,ij->k", trimers_formed, probs_above_threshold)
    trimers_below_threshold = np.einsum("ijk,ij->k", trimers_formed, 1-probs_above_threshold)
    ratio = trimers_above_threshold / (trimers_above_threshold + trimers_below_threshold + 1) # 加上1以避免0/0的错误
    p = (1 - (1 - self.prob_synapse)**trimers_above_threshold) * ratio # 一次接触能形成免疫突触的概率
    
    # 此阶段计算杀伤值
    _t = t.number(units.h)
    expected_synapses = _t * self.contact_freq * system.c[_effector, _compartments] * p # 期望免疫突触数量（对每个靶细胞来说），是长n_compartments的向量
    synapses = np.random.poisson(expected_synapses, size = self.hp.shape) # 实际免疫突触数量，是 n x n_compartments 的矩阵
    damages =  np.random.normal(synapses * self.damage_per_hit, np.sqrt(synapses) * self.cv_damage_per_hit)
    damages = np.maximum(0.0, damages)
    
    # 结算效果
    alives_prev = (self.hp > 0).mean(axis = 0)
    #assert alives_prev.min() > 0.1, "CytotoxicityProcess has dropped below 10 percent resolution!"
    self.hp += _t * self.repair - damages; self.hp = np.minimum(1.0, self.hp); self.hp[self.hp <= 0] = -np.inf
    alives = (self.hp > 0).mean(axis = 0)
    deaths = 1 - np.divide(alives, alives_prev, out = np.zeros_like(alives), where = alives_prev!=0) # 写成 1 - alives/alives_prev 会有0/0的情况发生
    system._cell_death(_affected, _compartments, deaths)
    '''
    if system.t > 11.9 and system.t < 12.1:
      print(f"12h杀伤监测: effective_trimers={trimers_above_threshold}\tp={p}\tdeaths={deaths}")
    '''
