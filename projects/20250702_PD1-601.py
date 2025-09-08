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


drug.set_affinity(",,", 1, "CD25", 3*units.nM, 4e-4/units.s)
drug.set_affinity(",,", 2, "CD132", 30*units.nM, 3e-3/units.s)


drug.set_trans("CD3", [",PSMA"], 5000/units.um**2, 9e-3/units.s)
drug.set_trans("PSMA", ["CD3,"], 500/units.um**2, 1e-4/units.s)
drug.set_internalization(["CD3,"], 0.02/units.h)
drug.set_internalization([",PSMA"], 0.1/units.h)
drug.set_synapse_probability(0.1)


