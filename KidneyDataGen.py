import argparse
from SaidmanCompatibleGenerator import SaidmanCompatibleGenerator
from SaidmanCompatibleGenerator import BJCSensitivityPool
from SaidmanCompatibleGenerator import BJCpool
from DistributionGenerator import *
import functions
import util
import json


K = 10
T = 10
filename = 'data.json'
pool = BJCSensitivityPool(T, K)
gen = DistributionGenerator()
matches = []

for i in range(K):
    matches.append([])
    for j in range(K):
        compatible_1 = functions.are_blood_compatible(pool.incompatiblePairs[i].bloodTypeDonor, pool.incompatiblePairs[j].bloodTypePatient) \
                  and (not pool.incompatiblePairs[i].saidman.isPositiveCrossmatch(pool.incompatiblePairs[j].patientCPRA))
        compatible_2 = functions.are_blood_compatible(pool.incompatiblePairs[j].bloodTypeDonor, pool.incompatiblePairs[i].bloodTypePatient) \
                    and (not pool.incompatiblePairs[j].saidman.isPositiveCrossmatch(pool.incompatiblePairs[i].patientCPRA))
        if not (compatible_1 and compatible_2): 
            matches[i].append(0)
            continue
        donor_rec_abo_comp = functions.are_blood_compatible(pool.incompatiblePairs[i].bloodTypeDonor, pool.incompatiblePairs[j].bloodTypePatient)
        donor_rec_HLA_B_mis = gen.gen_donor_rec_HLA_B_mis(0)
        donor_rec_HLA_DR_mis = gen.gen_donor_rec_HLA_DR_mis(0)
        donor_rec_weight_ratio = util.calculate_weight_ratio(pool.incompatiblePairs[i].donor_weight, pool.incompatiblePairs[j].rec_weight)

        LKDPI = util.calculate_lkdpi(pool.incompatiblePairs[i].donor_age, pool.incompatiblePairs[i].donor_afam, 
                                        pool.incompatiblePairs[i].donor_bmi, pool.incompatiblePairs[i].donor_cig_use,
                                          pool.incompatiblePairs[i].donor_sex, pool.incompatiblePairs[j].rec_sex,
                                          pool.incompatiblePairs[i].donor_sbp, donor_rec_abo_comp,
                                          0, pool.incompatiblePairs[i].donor_egfr, #assumed the donor and recip are NOT related
                                          donor_rec_HLA_B_mis, donor_rec_HLA_DR_mis,
                                          donor_rec_weight_ratio)
        matches[i].append(util.calculate_survival(LKDPI))

for i in range(T):
    matches.append([])
    for j in range(K):
        compatible_1 = functions.are_blood_compatible(pool.compatiblePairs[i].bloodTypeDonor, pool.incompatiblePairs[j].bloodTypePatient) \
                  and (not pool.compatiblePairs[i].saidman.isPositiveCrossmatch(pool.incompatiblePairs[j].patientCPRA))
        compatible_2 = functions.are_blood_compatible(pool.incompatiblePairs[j].bloodTypeDonor, pool.compatiblePairs[i].bloodTypePatient) \
                    and (not pool.incompatiblePairs[j].saidman.isPositiveCrossmatch(pool.compatiblePairs[i].patientCPRA))
        if not (compatible_1 and compatible_2): 
            matches[i+K].append(0)
            continue
        donor_rec_abo_comp = functions.are_blood_compatible(pool.incompatiblePairs[j].bloodTypeDonor, pool.compatiblePairs[i].bloodTypePatient)
        donor_rec_HLA_B_mis = gen.gen_donor_rec_HLA_B_mis(0)
        donor_rec_HLA_DR_mis = gen.gen_donor_rec_HLA_DR_mis(0)
        donor_rec_weight_ratio = util.calculate_weight_ratio(pool.incompatiblePairs[j].donor_weight, pool.compatiblePairs[i].rec_weight)

        LKDPI_IC = util.calculate_lkdpi(pool.incompatiblePairs[j].donor_age, pool.incompatiblePairs[j].donor_afam, 
                                        pool.incompatiblePairs[j].donor_bmi, pool.incompatiblePairs[j].donor_cig_use,
                                          pool.incompatiblePairs[j].donor_sex, pool.compatiblePairs[i].rec_sex,
                                          pool.incompatiblePairs[j].donor_sbp, donor_rec_abo_comp,
                                          0, pool.incompatiblePairs[j].donor_egfr, #assumed the donor and recip are NOT related
                                          donor_rec_HLA_B_mis, donor_rec_HLA_DR_mis,
                                          donor_rec_weight_ratio)
        
        if LKDPI_IC > pool.compatiblePairs[i].LKDPI:
            matches[i+K].append(0)
            matches[j].append(0)
            continue
        matches[j].append(util.calculate_survival(LKDPI_IC))
        donor_rec_abo_comp = functions.are_blood_compatible(pool.compatiblePairs[i].bloodTypeDonor, pool.incompatiblePairs[j].bloodTypePatient)
        donor_rec_HLA_B_mis = gen.gen_donor_rec_HLA_B_mis(0)
        donor_rec_HLA_DR_mis = gen.gen_donor_rec_HLA_DR_mis(0)
        donor_rec_weight_ratio = util.calculate_weight_ratio(pool.compatiblePairs[i].donor_weight, pool.incompatiblePairs[j].rec_weight)

        LKDPI_CI = util.calculate_lkdpi(pool.compatiblePairs[i].donor_age, pool.compatiblePairs[i].donor_afam, 
                                        pool.compatiblePairs[i].donor_bmi, pool.compatiblePairs[i].donor_cig_use,
                                          pool.compatiblePairs[i].donor_sex, pool.incompatiblePairs[j].rec_sex,
                                          pool.compatiblePairs[i].donor_sbp, donor_rec_abo_comp,
                                          0, pool.compatiblePairs[i].donor_egfr, #assumed the donor and recip are NOT related
                                          donor_rec_HLA_B_mis, donor_rec_HLA_DR_mis,
                                          donor_rec_weight_ratio)
        matches[i+K].append(util.calculate_survival(LKDPI_CI))
    for j in range(T):
        if j == i:
            matches[i+K].append(util.calculate_survival(pool.compatiblePairs[i].LKDPI))
        else:
            matches[i+K].append(0)


with open(filename, 'w') as f:
    f.write(json.dumps((K,T,matches)))


        
