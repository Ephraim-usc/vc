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
probody.set_affinity(",CD25,", 0, "PD1", 500/units.um**2, 1e-4/units.s)
probody.set_affinity(",,CD132", 0, "PD1", 500/units.um**2, 1e-4/units.s)

# 结合CD25后，对CD132亲和力增强100倍
probody.set_affinity(",CD25,", 2, "CD132", np.inf/units.um**2, 3e-3/units.s)
probody.set_affinity("PD1,CD25,", 2, "CD132", np.inf/units.um**2, 3e-3/units.s)
probody.set_affinity(",CD25,CD132", 0, "PD1", 500/units.um**2, 1e-4/units.s)

# 结合CD132后，对CD25亲和力增强10倍
probody.set_affinity(",,CD132", 1, "CD25", np.inf/units.um**2, 4e-4/units.s)
probody.set_affinity("PD1,,CD132", 1, "CD25", np.inf/units.um**2, 4e-4/units.s)


system = vc.InvivoSystem(plasma, lymph, [tumor, lung, liver, SI, spleen, other], [Tn, Teff], [probody, drug], [], halflives = [20*units.d, 20*units.d])
system.add_isodrug_transform(probody, drug, ["tumor", "plasma", "lymph", "SI", "lung", "other"], [0.094/units.d, 0.017/units.d, 0.017/units.d, 0.017/units.d, 0.017/units.d, 0.017/units.d])
system.add_drug_dose(drug, 3*units.mg/units.kg, molecular_weight = 100*units.kDa, body_weight = 70*units.kg)
system.run(24*7*units.h, t_step = 1/6 * units.h, verbose = True)
system.add_drug_dose(drug, 3*units.mg/units.kg, molecular_weight = 100*units.kDa, body_weight = 70*units.kg)
system.run(24*7*units.h, t_step = 1/6 * units.h, verbose = True)
system.add_drug_dose(drug, 3*units.mg/units.kg, molecular_weight = 100*units.kDa, body_weight = 70*units.kg)
system.run(24*7*units.h, t_step = 1/6 * units.h, verbose = True)
system.add_drug_dose(drug, 3*units.mg/units.kg, molecular_weight = 100*units.kDa, body_weight = 70*units.kg)
system.run(24*7*units.h, t_step = 1/6 * units.h, verbose = True)

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




