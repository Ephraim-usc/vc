import vc
from vc import units

from vc.cells import cell_lung, cell_SI, Tn, Teff, cell_PRAD, cell_spleen
from vc.human import plasma, lymph, tumor, lung, liver, SI, spleen, other
tumor.set_cell(cell_PRAD, density = 3E8/units.ml)

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import numpy as np
import pandas as pd


def print(system, filename, title):
  Teff_PD1_analytes = ['(Teff)PD1']
  Tn_PD1_analytes = ['(Tn)PD1']
  Teff_gamma_analytes = ['(Teff)unmasked[,,CD132]', '(Teff)unmasked[PD1,,CD132]', '(Teff)unmasked[PD1,CD25,CD132]']
  Tn_gamma_analytes = ['(Tn)unmasked[,,CD132]', '(Tn)unmasked[PD1,,CD132]', '(Tn)unmasked[PD1,CD25,CD132]']
  ts = [t/24 for t in system.history.data.keys()]
  area = Teff.area.number(units.um**2)
  
  fig, ax = plt.subplots()
  
  ax.plot(ts, system.history['y', Teff_gamma_analytes, 'tumor']["sum"]*area, label = "Tumor gamma bound", color = "tab:red")
  ax.plot(ts, system.history['y', Teff_PD1_analytes, 'tumor']["sum"]*area, label = "Tumor PD1 not-bound", color = "tab:red", linestyle = "dashed")
  
  ax.plot(ts, system.history['y', Tn_gamma_analytes, 'lung']["sum"]*area, label = "Lung gamma bound", color = "tab:blue")
  ax.plot(ts, system.history['y', Tn_PD1_analytes, 'lung']["sum"]*area, label = "Lung PD1 not-bound", color = "tab:blue", linestyle = "dashed")
  
  ax.set_xlabel("time (d)")
  ax.set_xlim(0, 28)
  ax.axvspan(0, 21, alpha=0.1, color = "grey")
  ax.axvspan(21, 42, alpha=0.2, color = "grey")
  ax.axvspan(42, 63, alpha=0.3, color = "grey")
  
  #ax.set_yscale("log")
  ax.set_ylabel("number")
  #ax.set_ylim(1e-6, 1e3)
  
  ax.set_title(title)
  ax.legend(loc = "upper right", prop={'size': 6})
  fig.tight_layout()
  plt.gcf().set_size_inches(8, 8)
  fig.savefig(filename, dpi = 300)
  plt.close(fig)












drug = vc.Drug("unmasked", 3)
drug.set_affinity(",,", 0, "PD1", 1*units.nM, 1e-4/units.s)
drug.set_affinity(",,", 1, "CD25", 3*units.nM, 4e-4/units.s)
drug.set_affinity(",,", 2, "CD132", 30*units.nM, 3e-3/units.s)
drug.set_internalization([",,", ",,CD132", ",CD25,", ",CD25,CD132", "PD1,,", "PD1,,CD132", "PD1,CD25,", "PD1,CD25,CD132",], 0.1/units.h)

# PD1对IL2R的CAB效应，按正常avidity计算
drug.set_affinity("PD1,,", 1, "CD25", 1500/units.um**2, 4e-4/units.s)
drug.set_affinity("PD1,,", 2, "CD132", 15000/units.um**2, 3e-3/units.s)

# IL2R对PD1的CAB效应，按正常avidity计算
drug.set_affinity(",CD25,", 0, "PD1", 500/units.um**2, 1e-4/units.s)
drug.set_affinity(",,CD132", 0, "PD1", 500/units.um**2, 1e-4/units.s)

# 结合CD25后，对CD132亲和力增强100倍
drug.set_affinity(",CD25,", 2, "CD132", 150/units.um**2, 3e-3/units.s)
drug.set_affinity("PD1,CD25,", 2, "CD132", 150/units.um**2, 3e-3/units.s)
drug.set_affinity(",CD25,CD132", 0, "PD1", 500/units.um**2, 1e-4/units.s)

# 结合CD132后，对CD25亲和力增强10倍
drug.set_affinity(",,CD132", 1, "CD25", 150/units.um**2, 4e-4/units.s)
drug.set_affinity("PD1,,CD132", 1, "CD25", 150/units.um**2, 4e-4/units.s)



probody = vc.Drug("masked", 3)
probody.set_affinity(",,", 0, "PD1", 1*units.nM, 1e-4/units.s)
probody.set_affinity(",,", 1, "CD25", np.inf*units.nM, 4e-4/units.s)
probody.set_affinity(",,", 2, "CD132", np.inf*units.nM, 3e-3/units.s)
probody.set_internalization([",,", ",,CD132", ",CD25,", ",CD25,CD132", "PD1,,", "PD1,,CD132", "PD1,CD25,", "PD1,CD25,CD132",], 0.1/units.h)

# PD1对IL2R的CAB效应，按正常avidity计算
probody.set_affinity("PD1,,", 1, "CD25", np.inf/units.um**2, 4e-4/units.s)
probody.set_affinity("PD1,,", 2, "CD132", np.inf/units.um**2, 3e-3/units.s)

# IL2R对PD1的CAB效应，按正常avidity计算
probody.set_affinity(",CD25,", 0, "PD1", 10/units.um**2, 1e-4/units.s)
probody.set_affinity(",,CD132", 0, "PD1", 10/units.um**2, 1e-4/units.s)

# 结合CD25后，对CD132亲和力增强100倍
probody.set_affinity(",CD25,", 2, "CD132", np.inf/units.um**2, 3e-3/units.s)
probody.set_affinity("PD1,CD25,", 2, "CD132", np.inf/units.um**2, 3e-3/units.s)
probody.set_affinity(",CD25,CD132", 0, "PD1", 10/units.um**2, 1e-4/units.s)

# 结合CD132后，对CD25亲和力增强10倍
probody.set_affinity(",,CD132", 1, "CD25", np.inf/units.um**2, 4e-4/units.s)
probody.set_affinity("PD1,,CD132", 1, "CD25", np.inf/units.um**2, 4e-4/units.s)















drug = vc.Drug("unmasked", 3)
drug.set_affinity(",,", 0, "PD1", 1*units.nM, 1e-4/units.s)
drug.set_affinity(",,", 1, "CD25", 3*units.nM, 4e-4/units.s)
drug.set_affinity(",,", 2, "CD132", 136*units.nM, 2e-2/units.s)
drug.set_internalization([",,", ",,CD132", ",CD25,", ",CD25,CD132", "PD1,,", "PD1,,CD132", "PD1,CD25,", "PD1,CD25,CD132",], 0.1/units.h)

# PD1对IL2R的CAB效应，按正常avidity计算
drug.set_affinity("PD1,,", 1, "CD25", 30/units.um**2, 4e-4/units.s)
drug.set_affinity("PD1,,", 2, "CD132", 1360/units.um**2, 2e-2/units.s)

# IL2R对PD1的CAB效应，按正常avidity计算
drug.set_affinity(",CD25,", 0, "PD1", 10/units.um**2, 1e-4/units.s)
drug.set_affinity(",,CD132", 0, "PD1", 10/units.um**2, 1e-4/units.s)

# 结合CD25后，对CD132亲和力增强100倍
drug.set_affinity(",CD25,", 2, "CD132", 13.6/units.um**2, 3e-3/units.s)
drug.set_affinity("PD1,CD25,", 2, "CD132", 13.6/units.um**2, 3e-3/units.s)
drug.set_affinity(",CD25,CD132", 0, "PD1", 10/units.um**2, 1e-4/units.s)

# 结合CD132后，对CD25亲和力增强10倍
drug.set_affinity(",,CD132", 1, "CD25", 3/units.um**2, 4e-4/units.s)
drug.set_affinity("PD1,,CD132", 1, "CD25", 3/units.um**2, 4e-4/units.s)


system = vc.InvivoSystem(plasma, lymph, [tumor, lung, liver, SI, spleen, other], [Tn, Teff], [drug], [], halflives = [20*units.d])
system.add_drug_dose(drug, 1*units.mg/units.kg, molecular_weight = 100*units.kDa, body_weight = 70*units.kg)
system.run(24*7*units.h, t_step = 1/6 * units.h, verbose = True)

Teff_gamma_analytes = ['(Teff)unmasked[,,CD132]', '(Teff)unmasked[PD1,,CD132]', '(Teff)unmasked[PD1,CD25,CD132]']
Tn_gamma_analytes = ['(Tn)unmasked[,,CD132]', '(Tn)unmasked[PD1,,CD132]', '(Tn)unmasked[PD1,CD25,CD132]']


system.history['y', Teff_gamma_analytes, 'tumor'].iloc[:,1:] * Teff.area.number(units.um**2)
system.history['y', Tn_gamma_analytes, 'tumor'].iloc[:,1:] * Tn.area.number(units.um**2)

system.history['y', Tn_gamma_analytes, 'lung'].iloc[:,1:] * Tn.area.number(units.um**2)









probody = vc.Drug("masked", 3)
probody.set_affinity(",,", 0, "PD1", 1*units.nM, 1e-4/units.s)
probody.set_affinity(",,", 1, "CD25", np.inf*units.nM, 4e-4/units.s)
probody.set_affinity(",,", 2, "CD132", np.inf*units.nM, 3e-3/units.s)
probody.set_internalization([",,", ",,CD132", ",CD25,", ",CD25,CD132", "PD1,,", "PD1,,CD132", "PD1,CD25,", "PD1,CD25,CD132",], 0.1/units.h)

# PD1对IL2R的CAB效应，按正常avidity计算
probody.set_affinity("PD1,,", 1, "CD25", np.inf/units.um**2, 4e-4/units.s)
probody.set_affinity("PD1,,", 2, "CD132", np.inf/units.um**2, 3e-3/units.s)

# IL2R对PD1的CAB效应，按正常avidity计算
probody.set_affinity(",CD25,", 0, "PD1", 10/units.um**2, 1e-4/units.s)
probody.set_affinity(",,CD132", 0, "PD1", 10/units.um**2, 1e-4/units.s)

# 结合CD25后，对CD132亲和力增强100倍
probody.set_affinity(",CD25,", 2, "CD132", np.inf/units.um**2, 3e-3/units.s)
probody.set_affinity("PD1,CD25,", 2, "CD132", np.inf/units.um**2, 3e-3/units.s)
probody.set_affinity(",CD25,CD132", 0, "PD1", 10/units.um**2, 1e-4/units.s)

# 结合CD132后，对CD25亲和力增强10倍
probody.set_affinity(",,CD132", 1, "CD25", np.inf/units.um**2, 4e-4/units.s)
probody.set_affinity("PD1,,CD132", 1, "CD25", np.inf/units.um**2, 4e-4/units.s)



system = vc.InvivoSystem(plasma, lymph, [tumor, lung, liver, SI, spleen, other], [Tn, Teff], [probody, drug], [], halflives = [20*units.d, 20*units.d])
system.add_isodrug_transform(probody, drug, ["tumor", "plasma", "lymph", "SI", "lung", "other"], [0.094/units.d, 0.017/units.d, 0.017/units.d, 0.017/units.d, 0.017/units.d, 0.017/units.d])
system.add_drug_dose(drug, 3*units.mg/units.kg, molecular_weight = 100*units.kDa, body_weight = 70*units.kg)
system.run(24*21*units.h, t_step = 2/6 * units.h, verbose = True)
system.add_drug_dose(drug, 3*units.mg/units.kg, molecular_weight = 100*units.kDa, body_weight = 70*units.kg)
system.run(24*21*units.h, t_step = 2/6 * units.h, verbose = True)
system.add_drug_dose(drug, 3*units.mg/units.kg, molecular_weight = 100*units.kDa, body_weight = 70*units.kg)
system.run(24*21*units.h, t_step = 2/6 * units.h, verbose = True)
print(system, "drug.png", "drug (IBI363-like)")


system = vc.InvivoSystem(plasma, lymph, [tumor, lung, liver, SI, spleen, other], [Tn, Teff], [probody, drug], [], halflives = [20*units.d, 20*units.d])
system.add_isodrug_transform(probody, drug, ["tumor", "plasma", "lymph", "SI", "lung", "other"], [0.094/units.d, 0.017/units.d, 0.017/units.d, 0.017/units.d, 0.017/units.d, 0.017/units.d])
system.add_drug_dose(probody, 3*units.mg/units.kg, molecular_weight = 100*units.kDa, body_weight = 70*units.kg)
system.run(24*21*units.h, t_step = 2/6 * units.h, verbose = True)
system.add_drug_dose(probody, 3*units.mg/units.kg, molecular_weight = 100*units.kDa, body_weight = 70*units.kg)
system.run(24*21*units.h, t_step = 2/6 * units.h, verbose = True)
system.add_drug_dose(probody, 3*units.mg/units.kg, molecular_weight = 100*units.kDa, body_weight = 70*units.kg)
system.run(24*21*units.h, t_step = 2/6 * units.h, verbose = True)
print(system, "probody.png", "probody (perfect masking)")


Teff_gamma_analytes = ['(Teff)unmasked[,,CD132]', '(Teff)unmasked[PD1,,CD132]', '(Teff)unmasked[PD1,CD25,CD132]']
Tn_gamma_analytes = ['(Tn)unmasked[,,CD132]', '(Tn)unmasked[PD1,,CD132]', '(Tn)unmasked[PD1,CD25,CD132]']

system.history['y', Teff_gamma_analytes, 'tumor'].iloc[:,1:] * Teff.area.number(units.um**2)
system.history['y', Tn_gamma_analytes, 'lung'].iloc[:,1:] * Tn.area.number(units.um**2)







system = vc.InvivoSystem(plasma, lymph, [tumor, lung, liver, SI, spleen, other], [Tn, Teff], [probody, drug], [], halflives = [20*units.d, 20*units.d])
system.add_isodrug_transform(probody, drug, ["tumor", "plasma", "lymph", "SI", "lung", "other"], [0.094/units.d, 0.017/units.d, 0.017/units.d, 0.017/units.d, 0.017/units.d, 0.017/units.d])

system.add_drug_dose(probody, 0.3*units.mg, molecular_weight = 100*units.kDa, body_weight = 70*units.kg)
system.run(24*7*units.h, t_step = 1/6 * units.h, verbose = True)
system.add_drug_dose(probody, 1.5*units.mg, molecular_weight = 100*units.kDa, body_weight = 70*units.kg)
system.run(24*7*units.h, t_step = 1/6 * units.h, verbose = True)
system.add_drug_dose(probody, 6*units.mg, molecular_weight = 100*units.kDa, body_weight = 70*units.kg)
system.run(24*7*units.h, t_step = 1/6 * units.h, verbose = True)
system.add_drug_dose(probody, 6*units.mg, molecular_weight = 100*units.kDa, body_weight = 70*units.kg)
system.run(24*7*units.h, t_step = 1/6 * units.h, verbose = True)




