import vc
from vc import units

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


from vc.cells import cell_lung, cell_SI, Tn, Teff, cell_PRAD, cell_spleen
from vc.human import plasma, lymph, tumor, lung, liver, SI, spleen, other
tumor.set_cell(cell_PRAD, density = 3E8/units.ml)

system = vc.InvivoSystem(plasma, lymph, [tumor, lung, liver, SI, spleen, other], [cell_lung, cell_SI, Tn, Teff, cell_PRAD, cell_spleen], [drug, probody], [], halflives = [70*units.h])
system.add_isodrug_transform(probody, drug, ["tumor", "plasma", "lymph", "SI", "lung", "other"], [0.3/units.d, 0.01/units.d, 0.02/units.d, 0.02/units.d, 0.02/units.d, 0.02/units.d])
system.add_drug_dose(probody, 0.1*units.mpk, molecular_weight = 100*units.kDa, body_weight = 70*units.kg)
system.run(24*units.h, t_step = 1/6 * units.h, verbose = True)

analytes = ['mTCE[,,]'] + system.get_all_analytes_of_drug_on_cell(probody, cell_PRAD) + ['TCE[,,]'] + system.get_all_analytes_of_drug_on_cell(drug, cell_PRAD)
system.history['x', analytes ,'tumor']
