############################# 导入库与定义单位 #############################
import math
import numpy as np
import pandas as pd
import itertools
import functools
import copy

from scipy.linalg import expm
from scipy.integrate import solve_ivp

from tqdm import tqdm
from time import time as tt

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

import unum
import unum.units as units

units.l = unum.new_unit('l', 1e-3 * units.m ** 3)
units.ml = unum.new_unit('ml', 1e-3 * units.l)
units.ul = unum.new_unit('ul', 1e-6 * units.l)
units.pl = unum.new_unit('pl', 1e-12 * units.l)

units.M = unum.new_unit('M', 1 * units.mol / units.l)
units.uM = unum.new_unit('uM', 1e-6 * units.mol / units.l)
units.nM = unum.new_unit('nM', 1e-9 * units.mol / units.l)

units.kDa = unum.new_unit('kDa', units.kg / units.mol)
units.avagadro = unum.new_unit('avagadro', 6.0221415e23 / units.mol)

units.mpk = units.mg / units.kg
units.upk = units.ug / units.kg

units.min = units.MIN

np.set_printoptions(suppress=True)


############################# 平方微分方程组求解器的包装 #############################

class RS: # linear and quadratic reaction system
  def __init__(self, n_analytes):
    self.active = False #体系内有反应则为active，否则运行时可以跳过
    self.refreshed = True #最后一次加入反应后如果refresh过则为refreshed，表明可以调用solve()
    self.n = n_analytes
    self.Q = np.zeros([n_analytes, n_analytes]) # linear term coefficients
    self.QQ = np.zeros([n_analytes, n_analytes, n_analytes]) # quadratic term coefficients
  
  def refresh(self):
    self.linear_i, self.linear_o = np.where(self.Q != 0)
    self.linear_k = self.Q[self.Q != 0]
    self.quadratic_i, self.quadratic_j, self.quadratic_o = np.where(self.QQ != 0)
    self.quadratic_k = self.QQ[self.QQ != 0]
    self.refreshed = True
  
  def add_reaction_forward(self, reactants, products, forward): # reactants can be of length 1 or 2, products can be of any length
    if len(reactants) == 1:
      self.Q[reactants, reactants] -= forward
      self.Q[reactants, products] += forward
    else:
      self.QQ[reactants[0], reactants[1], reactants] -= forward
      self.QQ[reactants[0], reactants[1], products] += forward
    self.active = True
    self.refreshed = False
  
  def add_reaction(self, reactants, products, forward, backward):
    self.active = True
    self.add_reaction_forward(reactants, products, forward)
    if backward > 0:
      self.add_reaction_forward(products, reactants, backward)
    self.active = True
    self.refreshed = False
  
  def rate(self, _, x):
    buffer = np.zeros(self.n)
    np.add.at(buffer, self.linear_o, self.linear_k * x[self.linear_i])
    np.add.at(buffer, self.quadratic_o, self.quadratic_k * x[self.quadratic_i] * x[self.quadratic_j])
    return buffer
  
  def jac(self, _, x):
    buffer = np.zeros((self.n, self.n))
    np.add.at(buffer, (self.linear_o, self.linear_i), self.linear_k)
    np.add.at(buffer, (self.quadratic_o, self.quadratic_i), self.quadratic_k * x[self.quadratic_j])
    np.add.at(buffer, (self.quadratic_o, self.quadratic_j), self.quadratic_k * x[self.quadratic_i])
    return buffer
  
  def solve(self, x, t):
    if not self.refreshed:
      self.refresh()
    return solve_ivp(self.rate, jac = self.jac, t_span = (0, t), y0 = x, t_eval=[t], method = "BDF", rtol = 1e-10, atol = 1e-8).y[:,0]


############################# 定义细胞和药物 #############################

# 一个细胞由面积（或半径）、表面蛋白及拷贝数、对杀伤的耐受参数等定义
class Cell:
  def __init__(self, name = None, diameter = None, area = None, markers = None, copies = None,
               damage_per_hit = 0.0, cv_damage_per_hit = 0.0, repair = 0.0/units.h):
    self.name = name
    if area is None:
      tmp = math.pi * diameter**2 / 4
      self.area = round(tmp.number(units.um**2), 4) * units.um**2
    else:
      self.area = area
    
    self.markers = markers
    self.copies = copies
    self.damage_per_hit = damage_per_hit
    self.cv_damage_per_hit = cv_damage_per_hit
    self.repair = repair


# 一个药物有n个结合域
# 通过set_affinity加入3D或in cis结合，从而维护一颗state树
# 通过set_trans加入细胞间形成trimer的反应
# 通过set_internalization设定各state下的内吞速率，将在初始化系统中被自动读入
class Drug: #以后要检查，不同药物之间对markers和solutes的定义有无矛盾！以及不同药物的name有无重复！
  def __init__(self, name, n = 0):
    self.name = name
    self.n = n
    self.markers = []
    self.solutes = []
    if n > 0:
      self.affinities = {","*(n-1): dict()}
    else:
      self.affinities = {}
    self.transes = {}
    self.internalization = {}
    self.synapse_probability = None
  
  def copy(self, newname):
    buffer = copy.deepcopy(self)
    buffer.name = newname
    return buffer
  
  def is_free(self, state):
    state_ = state.split(',')
    return all([ligand in self.solutes for ligand in state_ if ligand != ''])
  
  def is_bound(self, state): # 可以改成参照is_free的写法
    state_ = state.split(',')
    for ligand in state_:
      if ligand in self.markers:
        return True
    return False
  
  def set_affinity(self, state, binder, ligand, aff, off = None, soluble = False):
    # 检查和加入markers或solutes
    if soluble:
      assert ligand not in self.markers, f"{ligand} has been added as a marker, but now trying to add as a soluble!"
      if ligand not in self.solutes:
        self.solutes.append(ligand)
    else:
      assert ligand not in self.solutes, f"{ligand} has been added as a solute, but now trying to add as a marker!"
      if ligand not in self.markers:
        self.markers.append(ligand)
    
    # 检查单位是否正确
    if self.is_bound(state) and ligand in self.markers:
      assert aff.unit() == 1/units.um**2, "Cis affinity must be of unit 1/um**2!"
    else:
      assert aff.unit() == units.nM, "3D affinity must be of unit nM!"
    
    self.affinities[state][(binder, ligand)] = (aff, off/aff, off)
    state_ = state.split(','); state_[binder] = ligand; product = ','.join(state_)
    if product not in self.affinities.keys():
      self.affinities[product] = dict()
  
  # 药物各states结合对面细胞上的某个marker的能力
  def set_trans(self, marker, states, aff, off):
    if marker not in self.transes:
      self.transes[marker] = {}
    for state in states:
      self.transes[marker][state] = (aff, off/aff, off)
  
  # 药物的免疫突触质量
  def set_synapse_probability(self, synapse_probability):
    self.synapse_probability = synapse_probability
  
  def get_states_free(self, solutes = None):
    buffer = []
    for state in self.affinities.keys():
      state_ = state.split(',')
      is_included = True
      is_bound = False
      for ligand in state_:
        if ligand in self.markers:
          is_bound = True
        if solutes is not None:
          if ligand in self.solutes and ligand not in solutes: # 如果是solute，要求在solutes里
            is_included = False
      if is_included and not is_bound:
        buffer.append(state)
    return buffer
  
  def get_states_on_cell(self, cell, solutes = None):
    buffer = []
    for state in self.affinities.keys():
      state_ = state.split(',')
      is_included = True
      is_bound = False
      for ligand in state_:
        if ligand in self.markers and ligand not in cell.markers: # 如果是marker，要求在细胞的markers里
          is_included = False
        if solutes is not None:
          if ligand in self.solutes and ligand not in solutes: # 如果是solute，要求在solutes里
            is_included = False
        if ligand in self.markers: # 至少要有一个ligand是marker
          is_bound = True
      if is_included and is_bound:
        buffer.append(state)
    return buffer
  
  # 因为没有加上细胞名、药物名，所以不是关于实在的analytes，因此称为protoreaction
  # 只要反应物二者之一是游离的，那么就是3D反应，与产物是否游离无关
  # 如果反应物都不是游离的，那么就是2D反应，此时产物必定不是游离的
  def get_protoreactions(self):
    buffer_3d = []
    buffer_2d = []
    states_free = self.get_states_free()
    for state in self.affinities.keys():
      for binder, ligand in self.affinities[state].keys():
        state_ = state.split(','); state_[binder] = ligand; product = ','.join(state_)
        aff, on, off = self.affinities[state][(binder, ligand)]
        if state in states_free or ligand in self.solutes:
          buffer_3d.append((state, ligand, product, on, off))
        else:
          buffer_2d.append((state, ligand, product, on, off))
    return buffer_3d, buffer_2d
  
  # 在一个定义了细胞和solutes的体系下的所有reactions
  # 一个protoreaction如果在体系中找不到某ligand，则不会被实例化
  # 一个protoreaction如果在多个细胞上都能匹配（即能找到所有ligands），则会被多次实例化
  def get_reactions(self, cells, solutes):
    protoreactions_3d, protoreactions_2d = self.get_protoreactions()
    
    reactions_3d = []
    reactionses_2d = {cell:[] for cell in cells}
    for state, ligand, product, on, off in protoreactions_3d:
      if self.is_free(product): # 因为都是合成反应，所以只要产物是游离的，此反应就是游离的，就只需添加一次
        if not all([ligand in solutes for ligand in product.split(',') if ligand != '']): # 因为都是合成反应，所以只要产物的ligands都在体系内，此反应就存在于体系内
          continue
        analyte_state = f"{self.name}[{state}]"
        analyte_ligand = ligand
        analyte_product = f"{self.name}[{product}]"
        reactions_3d.append((analyte_state, analyte_ligand, analyte_product, on, off))
      else:
        for cell in cells:
          if not all([ligand in solutes + cell.markers for ligand in product.split(',') if ligand != '']): # 
            continue
          analyte_state = f"({cell.name}){self.name}[{state}]" if self.is_bound(state) else f"{self.name}[{state}]"
          analyte_ligand = f"({cell.name}){ligand}" if ligand in self.markers else ligand
          analyte_product = f"({cell.name}){self.name}[{product}]" if self.is_bound(product) else f"{self.name}[{product}]"
          reactions_3d.append((analyte_state, analyte_ligand, analyte_product, on, off))
    
    for state, ligand, product, on, off in protoreactions_2d:
      for cell in cells:
          if not all([ligand in solutes + cell.markers for ligand in product.split(',') if ligand != '']): # 
            continue
          analyte_state = f"({cell.name}){self.name}[{state}]" # 2D反应的反应物、ligand和产物肯定都在细胞上，就不用判断了
          analyte_ligand = f"({cell.name}){ligand}"
          analyte_product = f"({cell.name}){self.name}[{product}]"
          reactionses_2d[cell].append((analyte_state, analyte_ligand, analyte_product, on, off))
    
    return reactions_3d, reactionses_2d
  
  # 此方程用于在两个细胞间计算trans结合的系数
  # 简单来说，浓度乘积张量与此矩阵的乘积，即为trimer（在各compartment）的数量
  def _get_trans_mats(self, cell1, cell2, solutes):
    markers1 = cell1.markers
    markers2 = cell2.markers
    states1 = self.get_states_on_cell(cell1, solutes)
    states2 = self.get_states_on_cell(cell2, solutes)
    names1 = markers1 + [f"{self.name}[{state}]" for state in states1]
    names2 = markers2 + [f"{self.name}[{state}]" for state in states2]
    
    _mat_on = pd.DataFrame(0.0, index = names1, columns = names2)
    _mat_off = pd.DataFrame(0.0, index = names1, columns = names2)
    for marker in self.transes.keys():
      for state, value in self.transes[marker].items():
        aff, on, off = value
        if marker in markers1 and state in states2:
          _mat_on.loc[marker, f"{self.name}[{state}]"] = on.number(units.um**2/units.h)
          _mat_off.loc[marker, f"{self.name}[{state}]"] = off.number(1/units.h)
        if marker in markers2 and state in states1:
          _mat_on.loc[f"{self.name}[{state}]", marker] = on.number(units.um**2/units.h)
          _mat_off.loc[f"{self.name}[{state}]", marker] = off.number(1/units.h)
    
    idx1, idx2 = (_mat_on != 0).any(axis=1), (_mat_on != 0).any(axis=0)  # 以on为准，on为零的过程则实际上不发生
    _mat_on = _mat_on.loc[idx1, idx2]
    _mat_off = _mat_off.loc[idx1, idx2]
    names1 = _mat_on.index
    names2 = _mat_on.columns
    return names1, names2, _mat_on.values, _mat_off.values
  
  def set_internalization(self, states, rate):
    for state in states:
      self.internalization[state] = rate


############################# 定义系统 #############################

class History:
  def __init__(self, system, compartments, analytes, analyteses_bound, cells):
    self.system = system
    self.compartments = compartments
    self.analytes = analytes
    self.analyteses_bound = analyteses_bound
    self.analytes_bound = [analyte_bound for analytes_bound in analyteses_bound for analyte_bound in analytes_bound] # 把y叠成一个矩阵方便索引
    self.cells = cells
    self.data = dict()
  
  def save(self):
    t = self.system.t
    y = np.concatenate(self.system.y, axis=0)
    self.data[t] = {"x":self.system.x.copy(), "y":y, "c":self.system.c.copy()}
  
  def set(self, key, value):
    t = self.system.t
    self.data[t][key] = value
  
  def __getitem__(self, indices): # 维度分别为：key，analytes/analytes_bound/cells, compartments
    key, index, columns = indices
    if not isinstance(index, (list, tuple)):
      index = (index, )
    if not isinstance(columns, (list, tuple)):
      columns = (columns, )
    
    if key == "x":
      index_ = [self.analytes.index(_) for _ in index]
      columns_ = [self.compartments.index(_) for _ in columns]
    elif key == "y":
      index_ = [self.analytes_bound.index(_) for _ in index]
      columns_ = [self.compartments.index(_) for _ in columns]
    elif key == "c":
      index_ = [self.cells.index(_) for _ in index]
      columns_ = [self.compartments.index(_) for _ in columns]
    else:
      index_ = index
      columns_ = columns
    
    buffer = [[t, data[key][np.ix_(index_, columns_)].sum()] + data[key][np.ix_(index_, columns_)].ravel().tolist() for t, data in self.data.items()]
    if len(columns) == 1:
      names = ["t", "sum"] + [i for i in index]
    elif len(index) == 1:
      names = ["t", "sum"] + [j for j in columns]
    else:
      names = ["t", "sum"] + [str(i) + '|' + str(j) for j in columns for i in index]
    return pd.DataFrame(buffer, columns = names)




def decode_analyte(analyte):
  head = analyte.split(')')[0] + ')' if analyte[0] == '(' else ''
  tail = '[' + analyte.split('[')[-1] if analyte[-1] == ']' else ''
  substance = analyte[len(head):-len(tail)] if len(tail) else analyte[len(head):]
  return head, substance, tail


class System:
  def __init__(self, compartments, cells, drugs, solutes, volumes = None):
    self.compartments = compartments
    self.n_compartments = len(compartments)
    self.cells = cells; self.n_cells = len(cells)
    self.cell_areas = [cell.area.number(units.um**2) for cell in cells]
    self.drugs = drugs; self.n_drugs = len(drugs)
    self.solutes = solutes; self.n_solutes = len(solutes)
    if volumes is None: # 对没有flow的体系来说，体积一般不起作用，所以不是必需参数
      self.volumes = [1.0 for compartment in compartments]
    else:
      self.volumes = [volume.number(units.ml) for volume in volumes]
    self.register = None
    self.reports = dict() # 用来给各种process等活动提供报告输出空间
    
    analytes_free = solutes.copy()
    for drug in drugs:
      analytes_free += [f"{drug.name}[{state}]" for state in drug.get_states_free(solutes)]
    self.analytes_free = analytes_free
    self.n_analytes_free = len(analytes_free)
    
    analyteses_marker = [[] for cell in cells]
    for i, cell in enumerate(cells):
      analyteses_marker[i] += [f"({cell.name}){marker}" for marker in cell.markers]
    self.analyteses_marker = analyteses_marker
    self.ns_analyteses_marker = [len(analytes_marker) for analytes_marker in analyteses_marker]
    
    analyteses_bound = [[] for cell in cells]
    for i, cell in enumerate(cells):
      analyteses_bound[i] += [f"({cell.name}){marker}" for marker in cell.markers]
      for drug in drugs:
        analyteses_bound[i] += [f"({cell.name}){drug.name}[{state}]" for state in drug.get_states_on_cell(cell, solutes)]
    self.analyteses_bound = analyteses_bound
    self.ns_analyteses_bound = [len(analytes_bound) for analytes_bound in analyteses_bound]
    
    self.analytes = analytes_free + [analyte for analytes_bound in analyteses_bound for analyte in analytes_bound]
    self.analytes_free_ = [self.analytes.index(analyte_free) for analyte_free in analytes_free]
    self.analyteses_marker_ = [[self.analytes.index(analyte_marker) for analyte_marker in analytes_marker] for analytes_marker in analyteses_marker]
    self.analyteses_bound_ = [[self.analytes.index(analyte_bound) for analyte_bound in analytes_bound] for analytes_bound in analyteses_bound]
    self.n_analytes = len(self.analytes)
    
    # 初始化系统各方面动力学
    self.F = np.zeros([self.n_analytes, self.n_compartments, self.n_compartments], dtype = float) # flow matrix of free analytes between compartments, in 1/units.h （注意bound的analyte的flow是没有实际意义的，但作为功能可以保留）
    self.T = np.zeros([self.n_compartments, self.n_analytes, self.n_analytes], dtype = float) # transformation between analytes in each compartment, in 1/units.h
    self.M = np.zeros([self.n_cells, self.n_compartments, self.n_compartments], dtype = float) # migration matrix of cells between compartments, in 1/units.h
    self.RS_3d = [RS(self.n_analytes) for compartment in self.compartments] # 3D反应可能涉及所有analytes
    self.RS_2d = [[RS(len(self.analyteses_bound[i])) for i in range(self.n_cells)] for compartment in self.compartments] # 2D反应在各细胞上分别进行
    self.processes = []
    
    # 把drugs定义中内含的反应加入进来
    for drug in drugs:
      reactions_3d, reactionses_2d = drug.get_reactions(self.cells, self.solutes)
      for analyte_state, analyte_ligand, analyte_product, on, off in reactions_3d:
        self.add_universal_reaction_3d([analyte_state, analyte_ligand], [analyte_product], on, off)
      for cell, reactions_2d in reactionses_2d.items():
        for analyte_state, analyte_ligand, analyte_product, on, off in reactions_2d:
          self.add_universal_reaction_2d(cell, [analyte_state, analyte_ligand], [analyte_product], on, off)
    
    # 把drugs中定义的internalization加入进来
    for compartment in compartments:
      for drug in drugs:
        for analyte in self.analytes:
          head, substance, tail = decode_analyte(analyte)
          if head == '':
            continue
          if substance != drug.name:
            continue
          state = tail[1:-1]
          if state not in drug.internalization:
            continue
          rate = drug.internalization[state]
          ligands = [ligand for ligand in state.split(',') if ligand != '']
          products = [(ligand if ligand in self.solutes else head + ligand) for ligand in ligands]
          self.add_transform(compartment, analyte, products, rate)
    
    # 初始化各变量
    # compartments放在最后一个维度是为了方便各种向量运算
    self.t = 0
    self.x = np.zeros([self.n_analytes, self.n_compartments], dtype = float) # 3D concentration of analytes, in units.nM
    self.y = [np.zeros([len(analyteses_bound), self.n_compartments], dtype = float) for analyteses_bound in self.analyteses_bound] # 2D concentration of analytes, in 1/units.um**2
    self.c = np.zeros([self.n_cells, self.n_compartments], dtype = float) # concentration of cells, in 1/units.ml
    
    self.history = History(self, self.compartments, self.analytes, self.analyteses_bound, self.cells)
  
  # 下面两个函数方便找到drug在此系统中对应的所有analytes
  def get_all_analytes_of_drug(self, drug):
    return [analyte for analyte in self.analytes if analyte.split(')')[-1].split('[')[0] == drug.name]
  
  def get_all_analytes_of_drug_on_cell(self, drug, cell):
    return [f"({cell.name}){drug.name}[{state}]" for state in drug.get_states_on_cell(cell, self.solutes)]

  # 假定刚对x做了更新，现在要刷新y的数据
  # 1nM的抗原，若分布在1/ml的细胞上，细胞表面积为1um**2，则抗原二维密度为 1nM * ml * avagadro / um**2 = 602214150000/um**2
  def _refresh_y(self):
    for cell_, cell in enumerate(self.cells):
      # 原写法，会导致无细胞时出现0/0：self.y[cell_] = self.x[self.analyteses_bound_[cell_]] / self.c[cell_] / self.cell_areas[cell_] * 602214150000
      self.y[cell_] = np.divide(self.x[self.analyteses_bound_[cell_]] / self.cell_areas[cell_] * 602214150000, 
                                self.c[cell_], 
                                out = np.zeros_like(self.y[cell_]), 
                                where = self.c[cell_]!=0)
  
  # 假定刚对y做了更新，现在要刷新x的数据
  def _refresh_x(self):
    for cell_, cell in enumerate(self.cells):
      self.x[self.analyteses_bound_[cell_]] = self.y[cell_] * self.c[cell_] * self.cell_areas[cell_] * 1.66053886e-12
  
  def add_reaction_3d(self, compartment, reactants, products, forward, backward = None):
    compartment = self.compartments.index(compartment)
    reactants = [self.analytes.index(reactant) for reactant in reactants]
    products = [self.analytes.index(product) for product in products]
    forward = forward.number(units.nM**(1 - len(reactants)) / units.h)
    if backward is None:
      backward = 0.0
    else:
      backward = backward.number(units.nM**(1 - len(products)) / units.h)
    self.RS_3d[compartment].add_reaction(reactants, products, forward, backward)
  
  def add_reaction_2d(self, compartment, cell, reactants, products, forward, backward = None):
    compartment = self.compartments.index(compartment)
    cell = self.cells.index(cell)
    reactants = [self.analyteses_bound[cell].index(reactant) for reactant in reactants]
    products = [self.analyteses_bound[cell].index(product) for product in products]
    forward = forward.number((1/units.um**2)**(1 - len(reactants)) / units.h)
    if backward is None:
      backward = 0.0
    else:
      backward = backward.number((1/units.um**2)**(1 - len(products)) / units.h)
    self.RS_2d[compartment][cell].add_reaction(reactants, products, forward, backward)
  
  def add_universal_reaction_3d(self, reactants, products, forward, backward = None):
    for compartment in self.compartments:
      self.add_reaction_3d(compartment, reactants, products, forward, backward)
  
  def add_universal_reaction_2d(self, cell, reactants, products, forward, backward = None):
    for compartment in self.compartments:
      self.add_reaction_2d(compartment, cell, reactants, products, forward, backward)
  
  def add_flow(self, analyte, compartment_source, compartment_dest, rate):
    _compartment_source = self.compartments.index(compartment_source)
    _analyte = self.analytes.index(analyte)
    _rate = rate.number(1/units.h)
    self.F[_analyte, _compartment_source, _compartment_source] -= _rate
    if compartment_dest is not None:
      _compartment_dest = self.compartments.index(compartment_dest)
      self.F[_analyte, _compartment_source, _compartment_dest] += _rate * self.volumes[_compartment_source]/self.volumes[_compartment_dest]
  
  def add_drug_flow(self, drug, compartment_source, compartment_dest, rate):
    states_free = drug.get_states_free(self.solutes)
    analytes_free = [f"{drug.name}[{state}]" for state in states_free]
    for analyte_free in analytes_free:
      self.add_flow(analyte_free, compartment_source, compartment_dest, rate)
  
  def add_transform(self, compartment, analyte_source, analytes_dest, rate):
    _compartment = self.compartments.index(compartment)
    _rate = rate.number(1/units.h)
    _analyte_source = self.analytes.index(analyte_source)
    _analytes_dest = [self.analytes.index(analyte_dest) for analyte_dest in analytes_dest]
    self.T[_compartment, _analyte_source, _analyte_source] -= _rate
    for _analyte_dest in _analytes_dest:
      self.T[_compartment, _analyte_source, _analyte_dest] += _rate
  
  # isodrugs意思是结合的靶点完全一样，只是亲和力在数值上有差异，因此states之间有一一对应的关系
  # compartments和rates长度必须一致
  def add_isodrug_transform(self, drug1, drug2, compartments, rates):
    for analyte_source in self.analytes:
      if analyte_source.split(')')[-1].split('[')[0] != drug1.name: # 检验是否是属于drug1的analyte
        continue
      analyte_dest = analyte_source.replace(drug1.name + '[', drug2.name + '[')
      assert analyte_dest in self.analytes, f"Expected analyte_dest {analyte_dest} not found! Are these isodrugs?"
      for compartment, rate in zip(compartments, rates):
        self.add_transform(compartment, analyte_source, [analyte_dest], rate)
  
  def add_process(self, process):
    self.processes.append(process)
  
  def get_c(self, compartment, cell):
    compartment_ = self.compartments.index(compartment)
    cell_ = self.cells.index(cell)
    return self.c[cell_, compartment_] * 1/units.ml
  
  def add_cell(self, compartment, cell, value):
    compartment_ = self.compartments.index(compartment)
    cell_ = self.cells.index(cell)
    value_ = value.number(1/units.ml)
    
    # 如增加 1/ml 细胞，则增加 copy*1/ml/avagadro = copy*1.66053886e-12nM 的抗原
    self.x[self.analyteses_marker_[cell_], compartment_] += np.array(cell.copies) * value_ * 1.66053886e-12
    self.c[cell_, compartment_] += value_
    self._refresh_y()
  
  # deaths：长度为n_compartments的向量，表示各compartment的死亡比例
  def _cell_death(self, _cell, _compartments, deaths):
    alive = 1 - deaths
    self.c[_cell, _compartments] *= alive
    self.x[np.ix_(self.analyteses_bound_[_cell], _compartments)] *= alive
    self.y[_cell][:, _compartments] *= alive
  
  def get_x(self, compartment, analyte):
    compartment_ = self.compartments.index(compartment)
    analyte_ = self.analytes.index(analyte)
    return self.x[analyte_, compartment_] * units.nM
  
  def add_x(self, compartment, analyte, value):
    compartment_ = self.compartments.index(compartment)
    analyte_ = self.analytes.index(analyte)
    value_ = value.number(units.nM)
    self.x[analyte_, compartment_] += value_
    self._refresh_y()
  
  def add_drug(self, compartment, drug, value):
    compartment_ = self.compartments.index(compartment)
    analyte = f"{drug.name}[{drug.get_states_free(self.solutes)[0]}]"
    analyte_ = self.analytes.index(analyte)
    value_ = value.number(units.nM)
    self.x[analyte_, compartment_] += value_
  
  def _run_reactions(self, _t):
    for compartment_, RSs_2d in enumerate(self.RS_2d):
      for cell_, rs in enumerate(RSs_2d):
        if rs.active:
          self.y[cell_][:, compartment_] = self.RS_2d[compartment_][cell_].solve(self.y[cell_][:, compartment_], _t)
    self._refresh_x()
    
    for compartment_, rs in enumerate(self.RS_3d):
      if rs.active:
        self.x[:, compartment_] = self.RS_3d[compartment_].solve(self.x[:, compartment_], _t)
    self._refresh_y()
  
  # 将t作为t_delta进行一轮计算
  # 顺序是（原则上由快到慢）：2D reactions、3D reactions、transformations、flows、migrations
  def _run(self, _t):
    if self.register is not _t:
      _transforming_compartments = [_compartment for _compartment in range(self.n_compartments) if self.T[_compartment].any()]
      expmT = np.zeros([self.n_compartments, self.n_analytes, self.n_analytes], dtype = float)
      for _compartment in _transforming_compartments:
        expmT[_compartment] = expm(_t * self.T[_compartment])
      
      _flowing_analytes = [_analyte for _analyte in range(self.n_analytes) if self.F[_analyte].any()]
      expmF = np.zeros([self.n_analytes, self.n_compartments, self.n_compartments], dtype = float)
      for _analyte in _flowing_analytes:
        expmF[_analyte] = expm(_t * self.F[_analyte])
      
      self.register = _t
      self.load = _flowing_analytes, expmF, _transforming_compartments, expmT
    
    _flowing_analytes, expmF, _transforming_compartments, expmT = self.load
    
    self._run_reactions(_t)
    
    for _compartment in _transforming_compartments:
      self.x[:, _compartment] = np.dot(self.x[:, _compartment], expmT[_compartment])
    for _analyte in _flowing_analytes:
      self.x[_analyte] = np.dot(self.x[_analyte], expmF[_analyte])
    self._refresh_y() # 上面两步都不会使用到y，因此可以最后再更新y
    
    for process in self.processes:
      process(self, _t * units.h)
  
  def run(self, t, t_step = 1/60 * units.h, t_record = 1 * units.h, verbose = False):
    self.register = None
    _t = t.number(units.h)
    _t_start = self.t
    _t_end = _t_start + _t
    _t_step = t_step.number(units.h)
    _t_record = t_record.number(units.h)
    
    if verbose:
      pbar = tqdm(total = _t, unit = "h", bar_format = "{desc}: {percentage:3.0f}%|{bar}| {n:.2f}/{total_fmt} [{elapsed}<{remaining},  {rate_fmt}{postfix}]")
      pbar.update(0.0)
    while True:
      _t_prev = self.t
      _t_delta = min(_t_step, _t_end - self.t)
      
      self._run(_t_delta)
      self.t = min(self.t + _t_step, _t_end) #采用这种写法，而不是self.t + _t_delta，是为了避免加法的偏差，尽量保持一个“整数”
      
      if math.floor(self.t / _t_record) > math.floor(_t_prev / _t_record):
        self.history.save()
      if verbose:
        pbar.update(_t_delta)
      if math.isclose(self.t, _t_end, rel_tol = 0, abs_tol = 1e-9):
        break
    if verbose:
      pbar.close()


############################# 定义器官与InVivo系统 #############################

# 一个器官由体积和体液进出速度定义
# 器官中可初始化所含的细胞和溶质
class Organ:
  def __init__(self, name, volume):
    self.name = name
    self.volume = volume
    self.entry = 0.0/units.h
    self.exit_plasma = 0.0/units.h
    self.exit_lymph = 0.0/units.h
    self.cells = {}
    self.solutes = []
  
  def set_entry(self, rate):
    self.entry = rate
  
  def set_exit_plasma(self, rate):
    self.exit_plasma = rate
  
  def set_exit_lymph(self, rate):
    self.exit_lymph = rate
  
  def add_convection(self, plasma_flow = None, lymph_flow_ratio = None, plasma_reflection = None, lymph_reflection = None, ABC = None, entering_time = None):
    if entering_time is not None:
      self.entry += (math.log(2)/entering_time).cast_unit(1/units.h)
      self.exit_lymph += (math.log(2)/entering_time / ABC).cast_unit(1/units.h)
      return
    if plasma_reflection is None:
      plasma_reflection = 1 - ABC * (1 - lymph_reflection)  # (1 - plasma_reflection) / (1 - lymph_reflection) == ABC
    self.entry += (plasma_flow * lymph_flow_ratio * (1 - plasma_reflection) / self.volume).cast_unit(1/units.h)
    self.exit_lymph += (plasma_flow * lymph_flow_ratio * (1 - lymph_reflection) / self.volume).cast_unit(1/units.h)
  
  def add_diffusion(self, plasma_interstitial_ratio, capillary_radius, capillary_permeability):
    rate = 2 * plasma_interstitial_ratio * capillary_permeability / capillary_radius
    self.entry += rate.cast_unit(1/units.h)
    self.exit_plasma += rate.cast_unit(1/units.h)
  
  def set_cell(self, cell, total = None, density = None):
    if density is None:
      density = total / self.volume
    self.cells[cell] = density
  
  def get_ABC(self):
    return (self.entry / (self.exit_plasma + self.exit_lymph)).number()


# InVivo系统是系统的子类，要求必须有血液和淋巴
# 作为系统的子类，需要初始化包含的细胞、药物和溶质，并且额外的可以初始化药物的半衰期
# 根据各器官定义的体液流动速度自动建立flow
# 根据各器官定义的初始化自动添加细胞和溶质
class InvivoSystem(System):
  def __init__(self, plasma, lymph, organs, cells, drugs, solutes, halflives = None):
    compartments = [plasma.name, lymph.name] + [organ.name for organ in organs]
    volumes = [plasma.volume, lymph.volume] + [organ.volume for organ in organs]
    self.plasma = plasma
    self.lymph = lymph
    self.organs = organs
    System.__init__(self, compartments, cells, drugs, solutes, volumes)
    
    # 加入各器官的流动
    for organ in self.organs:
      for analyte in self.analytes_free:
        self.add_flow(analyte, self.plasma.name, organ.name, organ.entry * organ.volume / self.plasma.volume)
        self.add_flow(analyte, organ.name, self.plasma.name, organ.exit_plasma)
        self.add_flow(analyte, organ.name, self.lymph.name, organ.exit_lymph)
        self.add_flow(analyte, self.lymph.name, self.plasma.name, organ.exit_lymph * organ.volume / self.lymph.volume)
    
    # 加入各药物的清除
    if halflives:
      for drug, halflife in zip(drugs, halflives):
        rate = math.log(2)/halflife # exp(- rate * halflife) = 1/2
        for compartment in self.compartments:
          self.add_drug_flow(drug, compartment, None, rate)
    
    # 加入各器官中的细胞
    for organ in [plasma, lymph] + organs:
      for cell, density in organ.cells.items():
        if cell not in self.cells:
          continue
        self.add_cell(organ.name, cell, density)
  
  # 接受mpk或mg两种单位的剂量打入
  def add_drug_dose(self, drug, value, molecular_weight, body_weight): # 标准体重：小鼠28g，猴6.2kg，人71kg
    if value.unit() == units.mg:
      conc = value / molecular_weight / self.plasma.volume
    elif value.unit() == units.mpk:
      conc = value * body_weight / molecular_weight / self.plasma.volume
    else:
      print("unrecognized dose unit!")
      conc = 0*units.nM
    self.add_drug(self.plasma.name, drug, conc)
