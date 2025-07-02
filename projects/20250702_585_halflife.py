import vc
from vc import units

drug = vc.Drug("TCE", 4)
drug.set_affinity(",,,", 0, "CD3", 10*units.nM, 9e-3/units.s)
drug.set_affinity(",,,", 1, "PD1", 1*units.nM, 1e-4/units.s)
drug.set_affinity(",,,", 2, "PSMA", 1*units.nM, 1e-4/units.s)
drug.set_affinity(",,,", 3, "STEAP1", 10*units.nM, 1e-3/units.s)
drug.set_affinity(",PD1,,", 0, "CD3", 5000/units.um**2, 9e-3/units.s)
drug.set_affinity("CD3,,,", 1, "PD1", 500/units.um**2, 1e-4/units.s)
drug.set_affinity(",,,STEAP1", 2, "PSMA", 500/units.um**2, 1e-4/units.s)
drug.set_affinity(",,PSMA,", 3, "STEAP1", 5000/units.um**2, 1e-3/units.s)
drug.set_trans("CD3", [",,PSMA,", ",,,STEAP1", ",,PSMA,STEAP1"], 5000/units.um**2, 9e-3/units.s)
drug.set_trans("PSMA", ["CD3,,,", "CD3,PD1,,"], 500/units.um**2, 1e-4/units.s)
drug.set_trans("STEAP1", ["CD3,,,", "CD3,PD1,,"], 5000/units.um**2, 1e-3/units.s)
drug.set_internalization(["CD3,,,"], 0.02/units.h)
drug.set_internalization([",PD1,,"], 0.02/units.h)
drug.set_internalization(["CD3,PD1,,"], 0.02/units.h)
drug.set_internalization([",,PSMA,"], 0.1/units.h)
drug.set_internalization([",,,STEAP1"], 0.1/units.h)
drug.set_internalization([",,PSMA,STEAP1"], 0.1/units.h)
drug.set_synapse_probability(0.1)

probody = vc.Drug("mTCE", 4)
probody.set_affinity(",,,", 0, "CD3", 200*units.nM, 9e-3/units.s)
probody.set_affinity(",,,", 1, "PD1", 1*units.nM, 1e-4/units.s)
probody.set_affinity(",,,", 2, "PSMA", 20*units.nM, 1e-4/units.s)
probody.set_affinity(",,,", 3, "STEAP1", 10*units.nM, 1e-3/units.s)
probody.set_affinity(",PD1,,", 0, "CD3", 100000/units.um**2, 9e-3/units.s)
probody.set_affinity("CD3,,,", 1, "PD1", 500/units.um**2, 1e-4/units.s)
probody.set_affinity(",,,STEAP1", 2, "PSMA", 10000/units.um**2, 1e-4/units.s)
probody.set_affinity(",,PSMA,", 3, "STEAP1", 5000/units.um**2, 1e-3/units.s)
probody.set_trans("CD3", [",,PSMA,", ",,,STEAP1", ",,PSMA,STEAP1"], 10000/units.um**2, 9e-3/units.s)
probody.set_trans("PSMA", ["CD3,,,", "CD3,PD1,,"], 10000/units.um**2, 1e-4/units.s)
probody.set_trans("STEAP1", ["CD3,,,", "CD3,PD1,,"], 5000/units.um**2, 1e-3/units.s)
probody.set_internalization(["CD3,,,"], 0.02/units.h)
probody.set_internalization([",PD1,,"], 0.02/units.h)
probody.set_internalization(["CD3,PD1,,"], 0.02/units.h)
probody.set_internalization([",,PSMA,"], 0.1/units.h)
probody.set_internalization([",,,STEAP1"], 0.1/units.h)
probody.set_internalization([",,PSMA,STEAP1"], 0.1/units.h)
probody.set_synapse_probability(0.1)



for psma in [3210330]:
  for aff in [0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100]:
    cell_PRAD.copies =[psma, 198752]
    
    drug.set_affinity(",,,", 1, "PD1", aff*units.nM, aff*1e-4/units.s)
    drug.set_affinity("CD3,,,", 1, "PD1", aff*500/units.um**2, aff*1e-4/units.s)
    
    probody.set_affinity(",,,", 1, "PD1", aff*units.nM, aff*1e-4/units.s)
    probody.set_affinity("CD3,,,", 1, "PD1", aff*500/units.um**2, aff*1e-4/units.s)
    
    system = InvivoSystem(plasma, lymph, [tumor, lung, liver, SI, spleen, other], [cell_lung, cell_SI, Tn, Teff, cell_PRAD, cell_spleen], [drug, probody], [], halflives = [70*units.h])
    system.add_isodrug_transform(probody, drug, ["tumor", "plasma", "lymph", "SI", "lung", "other"], [0.3/units.d, 0.01/units.d, 0.02/units.d, 0.02/units.d, 0.02/units.d, 0.02/units.d])
    system.add_drug_dose(probody, 0.1*units.mpk, molecular_weight = 100*units.kDa, body_weight = 70*units.kg)
    p_Tn_tumor = ComputeTrimersProcess(system, Tn, cell_PRAD, drug, ["tumor"])
    p_Tn_SI = ComputeTrimersProcess(system, Tn, cell_SI, drug, ["SI"])
    p_Tn_lung = ComputeTrimersProcess(system, Tn, cell_lung, drug, ["lung"])
    p_Teff_tumor = ComputeTrimersProcess(system, Teff, cell_PRAD, drug, ["tumor"])
    for p in p_Tn_tumor, p_Tn_SI, p_Tn_lung, p_Teff_tumor:
      system.add_process(p)
    system.run(12*units.h, t_step = 1/6 * units.h, verbose = True)
    
    #print(system.get_total_drug_in_compartment(["tumor", "plasma", "lymph", "SI", "lung", "spleen"], probody))
    #print(system.get_total_drug_in_compartment(["tumor", "plasma", "lymph", "SI", "lung", "spleen"], drug))
    print(system.get_total_drug_in_compartment(["tumor", "plasma", "lymph", "SI", "lung", "spleen"], probody) + system.get_total_drug_in_compartment(["tumor", "plasma", "lymph", "SI", "lung", "spleen"], drug))


