import vc
from vc import units

from vc.cells import cell_lung, cell_SI, Tn, Teff, cell_PRAD, cell_spleen
from vc.human import plasma, lymph, tumor, lung, liver, SI, spleen, other
tumor.set_cell(cell_PRAD, density = 3E8/units.ml)



drug = vc.Drug("TCE", 3)
drug.set_affinity(",,", 0, "CD3", 10*units.nM, 9e-3/units.s)
drug.set_affinity(",,", 1, "PSMA", 1*units.nM, 1e-4/units.s)
drug.set_affinity(",,", 2, "STEAP1", 10*units.nM, 1e-3/units.s)
drug.set_affinity(",,STEAP1", 1, "PSMA", 500/units.um**2, 1e-4/units.s)
drug.set_affinity(",PSMA,", 2, "STEAP1", 5000/units.um**2, 1e-3/units.s)
drug.set_trans("CD3", [",PSMA,", ",,STEAP1", ",PSMA,STEAP1"], 5000/units.um**2, 9e-3/units.s)
drug.set_trans("PSMA", ["CD3,,"], 500/units.um**2, 1e-4/units.s)
drug.set_trans("STEAP1", ["CD3,,"], 5000/units.um**2, 1e-3/units.s)
drug.set_internalization(["CD3,,"], 0.02/units.h)
drug.set_internalization([",PSMA,"], 0.1/units.h)
drug.set_internalization([",,STEAP1"], 0.1/units.h)
drug.set_internalization([",PSMA,STEAP1"], 0.1/units.h)
drug.set_synapse_probability(0.1)

probody = vc.Drug("mTCE", 3)
probody.set_affinity(",,", 0, "CD3", 200*units.nM, 9e-3/units.s)
probody.set_affinity(",,", 1, "PSMA", 20*units.nM, 1e-4/units.s)
probody.set_affinity(",,", 2, "STEAP1", 10*units.nM, 1e-3/units.s)
probody.set_affinity(",,STEAP1", 1, "PSMA", 10000/units.um**2, 1e-4/units.s)
probody.set_affinity(",PSMA,", 2, "STEAP1", 5000/units.um**2, 1e-3/units.s)
probody.set_trans("CD3", [",PSMA,", ",,STEAP1", ",PSMA,STEAP1"], 10000/units.um**2, 9e-3/units.s)
probody.set_trans("PSMA", ["CD3,,"], 10000/units.um**2, 1e-4/units.s)
probody.set_trans("STEAP1", ["CD3,,"], 5000/units.um**2, 1e-3/units.s)
probody.set_internalization(["CD3,,"], 0.02/units.h)
probody.set_internalization([",PSMA,"], 0.1/units.h)
probody.set_internalization([",,STEAP1"], 0.1/units.h)
probody.set_internalization([",PSMA,STEAP1"], 0.1/units.h)
probody.set_synapse_probability(0.1)


import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import numpy as np
import pandas as pd

def print_pk(system, filename):
  analytes_probody = system.get_all_analytes_of_drug(probody)
  analytes_drug = system.get_all_analytes_of_drug(drug)
  analytes = analytes_probody + analytes_drug
  ts = system.history.data.keys()
  
  fig, ax = plt.subplots()
  
  ax.plot(ts, system.history['x', analytes, 'plasma']["sum"], label = "plasma total drug", color = "tab:red")
  ax.plot(ts, system.history['x', analytes, 'tumor']["sum"], label = "tumor total drug", color = "tab:blue")
  
  ax.plot(ts, system.history['x', analytes_drug, 'plasma']["sum"], label = "plasma uTCE", color = "tab:red", linestyle = "dashed")
  ax.plot(ts, system.history['x', analytes_drug, 'tumor']["sum"], label = "tumor uTCE", color = "tab:blue", linestyle = "dashed")
  
  ax.set_yscale("log")
  
  ax.legend(loc = "upper right", prop={'size': 6})
  fig.tight_layout()
  plt.gcf().set_size_inches(8, 7)
  fig.savefig("tmp.png", dpi = 300)
  plt.close(fig)






drug = vc.Drug("TCE", 2)
drug.set_affinity(",", 0, "CD3", 10*units.nM, 9e-3/units.s)
drug.set_affinity(",", 1, "PSMA", 1*units.nM, 1e-4/units.s)
drug.set_trans("CD3", [",PSMA"], 5000/units.um**2, 9e-3/units.s)
drug.set_trans("PSMA", ["CD3,"], 500/units.um**2, 1e-4/units.s)
drug.set_internalization(["CD3,"], 0.02/units.h)
drug.set_internalization([",PSMA"], 0.1/units.h)
drug.set_synapse_probability(0.1)

probody = vc.Drug("mTCE", 2)
probody.set_affinity(",", 0, "CD3", 200*units.nM, 9e-3/units.s)
probody.set_affinity(",", 1, "PSMA", 20*units.nM, 1e-4/units.s)
probody.set_trans("CD3", [",PSMA"], 10000/units.um**2, 9e-3/units.s)
probody.set_trans("PSMA", ["CD3,"], 10000/units.um**2, 1e-4/units.s)
probody.set_internalization(["CD3,"], 0.02/units.h)
probody.set_internalization([",PSMA"], 0.1/units.h)
probody.set_synapse_probability(0.1)

system = vc.InvivoSystem(plasma, lymph, [tumor, lung, liver, SI, spleen, other], [cell_lung, cell_SI, Tn, Teff, cell_PRAD, cell_spleen], [probody, drug], [], halflives = [70*units.h, 1*units.h])
system.add_isodrug_transform(probody, drug, ["tumor", "plasma", "lymph", "SI", "lung", "other"], [0.3/units.d, 0.01/units.d, 0.02/units.d, 0.02/units.d, 0.02/units.d, 0.02/units.d])
system.add_drug_dose(probody, 0.1*units.mpk, molecular_weight = 100*units.kDa, body_weight = 70*units.kg)
system.run(24*14*units.h, t_step = 1/6 * units.h, verbose = True)

print_pk(system, "JANX007.png")

