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


drug = vc.Drug("drug", 3)
drug.set_affinity(",,", 0, "PD1", 1*units.nM, 1e-4/units.s)
drug.set_affinity(",,", 1, "CD25", 3*units.nM, 4e-4/units.s)
drug.set_affinity(",,", 2, "CD132", 30*units.nM, 3e-3/units.s)

# PD1对IL2R的CAB效应，按正常avidity计算
drug.set_affinity("PD1,,", 1, "CD25", 1500/units.mm**2, 4e-4/units.s)
drug.set_affinity("PD1,,", 2, "CD132", 15000/units.mm**2, 3e-3/units.s)

# IL2R对PD1的CAB效应，按正常avidity计算
drug.set_affinity(",CD25,", 0, "PD1", 500/units.mm**2, 1e-4/units.s)
drug.set_affinity(",,CD132", 0, "PD1", 500/units.mm**2, 1e-4/units.s)
drug.set_affinity(",CD25,CD132", 0, "PD1", 500/units.mm**2, 1e-4/units.s)

# 结合CD25后，对CD132亲和力增强100倍
drug.set_affinity(",CD25,", 2, "CD132", 150/units.mm**2, 3e-3/units.s)
drug.set_affinity("PD1,CD25,", 2, "CD132", 150/units.mm**2, 3e-3/units.s)

# 结合CD132后，对CD25亲和力增强10倍
drug.set_affinity(",,CD132", 1, "CD25", 150/units.mm**2, 4e-4/units.s)
drug.set_affinity("PD1,,CD132", 1, "CD25", 150/units.mm**2, 4e-4/units.s)

drug.set_internalization([",,", ",,CD132", ",CD25,", ",CD25,CD132", "PD1,,", "PD1,,CD132", "PD1,CD25,", "PD1,CD25,CD132",], 0.1/units.h)
drug.set_synapse_probability(0.1)


