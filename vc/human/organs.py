from vc import *
from vc.cells import *

plasma = Organ("plasma", 3126 * units.ml)
plasma.set_cell(Tn, total = 7.9E+09)

lymph = Organ("lymph", 3126 * units.ml)
lymph.set_cell(Tn, total = 3.6E+11)

lung = Organ("lung", 300 * units.ml)
lung.add_convection(plasma_flow = 181913 * units.ml/units.h, lymph_flow_ratio = 0.002, plasma_reflection = 0.9, lymph_reflection = 0.2)
lung.set_cell(cell_lung, total = 2.36E+11 * 0.5)
lung.set_cell(Tn, total = 1.3E+10)

liver = Organ("liver", 429 * units.ml)
liver.add_convection(plasma_flow = 13210 * units.ml/units.h, lymph_flow_ratio = 0.002, plasma_reflection = 0.85, lymph_reflection = 0.2)
liver.set_cell(Tn, total = 7.9E+09)

SI = Organ("SI", 67.1 * units.ml)
SI.add_convection(plasma_flow = 12368 * units.ml/units.h, lymph_flow_ratio = 0.002, plasma_reflection = 0.9, lymph_reflection = 0.2)
SI.set_cell(cell_SI, total = 7.2e11 * 0.5)
SI.set_cell(Tn, total = 1.8E+10)

spleen = Organ("spleen", 44.3 * units.ml)
spleen.add_convection(plasma_flow = 6343 * units.ml/units.h, lymph_flow_ratio = 0.002, plasma_reflection = 0.85, lymph_reflection = 0.2)
spleen.set_cell(cell_spleen, total = 1.83E+12 * 0.5)
spleen.set_cell(Tn, total = 221*3e8)

other = Organ("other", 5000 * units.ml)
other.add_convection(plasma_flow = 100000 * units.ml/units.h, lymph_flow_ratio = 0.002, plasma_reflection = 0.95, lymph_reflection = 0.2)
other.set_cell(Tn, total = 1.1E+10)

tumor = Organ("tumor", 170 * units.ul * 0.55)
tumor.add_convection(plasma_flow = 12.7/units.h * 170*units.ul, lymph_flow_ratio = 0.002, plasma_reflection = 0.842, lymph_reflection = 0.2)
tumor.add_diffusion(0.07/0.55, capillary_radius = 10 * units.um, capillary_permeability = 3e-7 * units.cm/units.s)
tumor.set_cell(Tn, density = 3E8/units.ml * 0.1)
tumor.set_cell(Teff, density = 3E8/units.ml * 0.1)
#tumor.set_cell(cell_PRAD, density = 3E8/units.ml)
