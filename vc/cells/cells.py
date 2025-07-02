import vc

# 免疫细胞
Tn = vc.Cell("Tn", diameter = 6*units.um, markers = ["CD3", "PD1"], copies = [50000, 10000])
Teff = vc.Cell("Teff", diameter = 6*units.um, markers = ["CD3", "PD1"], copies = [50000, 50000])

# 组织细胞
cell_lung = vc.Cell("lung", diameter = 50*units.um, markers = ["PSMA", "STEAP1"], copies = [22069, 39737])
cell_SI = vc.Cell("SI", diameter = 20*units.um, markers = ["PSMA", "STEAP1"], copies = [255438, 2802])
cell_spleen = vc.Cell("spleen", diameter = 20*units.um, markers = ["PSMA", "STEAP1"], copies = [0, 0])

# 细胞株
cell_HT29 = vc.Cell("HT29", diameter = 18*units.um, markers = ["EGFR", "CA9"], copies = [50000, 50000], damage_per_hit = 0.33, cv_damage_per_hit = 0.33, repair = 0.05/units.h)

# 适应症典型细胞
cell_PRAD = vc.Cell("PRAD", diameter = 15*units.um, markers = ["PSMA", "STEAP1"], copies = [3210330, 198752])
