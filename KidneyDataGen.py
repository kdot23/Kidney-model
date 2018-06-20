import argparse
from SaidmanCompatibleGenerator import SaidmanCompatibleGenerator
from SaidmanCompatibleGenerator import BJCSensitivityPool
from SaidmanCompatibleGenerator import BJCpool
from DistributionGenerator import *
import functions
import util
import json

parser = argparse.ArgumentParser(description="Generates Donor recipient pairs and a quality pool for optimization")
parser.add_argument('-K', '--num_incompatible', default=100, dest='K', type = int)
parser.add_argument('-T', '--num_compatible', default=100, dest='T', type = int)
parser.add_argument('-o', '--output', default='data.json')

args = parser.parse_args()

K = args.K
T = args.T
filename = args.output
pool = BJCSensitivityPool(T, K)
gen = DistributionGenerator()
matches = []
demo = []

for i in range(T):
    matches.append([])
    matches[i].append(util.calculate_survival(pool.compatiblePairs[i].LKDPI))
    demo.append([])
    demo[i] = (int(pool.compatiblePairs[i].bloodTypePatient == 0), int( pool.compatiblePairs[i].bloodTypePatient == 1 ), \
        int( pool.compatiblePairs[i].bloodTypePatient == 2 ), int( pool.compatiblePairs[i].bloodTypePatient == 3 ), int( pool.compatiblePairs[i].bloodTypeDonor == 0 ), \
        int( pool.compatiblePairs[i].bloodTypeDonor == 1 ), int( pool.compatiblePairs[i].bloodTypeDonor == 2 ), int( pool.compatiblePairs[i].bloodTypeDonor == 3 ), \
        pool.compatiblePairs[i].donor_afam, pool.compatiblePairs[i].donor_age, pool.compatiblePairs[i].donor_sex[0], \
        pool.compatiblePairs[i].donor_cig_use[0], pool.compatiblePairs[i].rec_sex[0], pool.compatiblePairs[i].donor_weight, pool.compatiblePairs[i].rec_weight, \
        pool.compatiblePairs[i].donor_bmi)
    for j in range(K):
        compatible_1 = functions.are_blood_compatible(pool.compatiblePairs[i].bloodTypeDonor, pool.incompatiblePairs[j].bloodTypePatient) \
                  and (not pool.compatiblePairs[i].saidman.isPositiveCrossmatch(pool.incompatiblePairs[j].patientCPRA))
        compatible_2 = functions.are_blood_compatible(pool.incompatiblePairs[j].bloodTypeDonor, pool.compatiblePairs[i].bloodTypePatient) \
                    and (not pool.incompatiblePairs[j].saidman.isPositiveCrossmatch(pool.compatiblePairs[i].patientCPRA))
        if not (compatible_1 and compatible_2): 
            matches[i].append(0)
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
            matches[i].append(0)
            continue
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
        matches[i].append(util.calculate_survival(LKDPI_CI)+util.calculate_survival(LKDPI_IC))


for i in range(K):
    matches.append([])
    matches[i+T].append(0)
    demo.append([])
    demo[i+T] = (int( pool.incompatiblePairs[i].bloodTypePatient == 0 ), int( pool.incompatiblePairs[i].bloodTypePatient == 1 ), int( pool.incompatiblePairs[i].bloodTypePatient == 2 ), \
        int( pool.incompatiblePairs[i].bloodTypePatient == 3 ), int( pool.incompatiblePairs[i].bloodTypeDonor == 0 ), int( pool.incompatiblePairs[i].bloodTypeDonor == 1 ), pool.incompatiblePairs[i].bloodTypeDonor == 2, \
        int( pool.incompatiblePairs[i].bloodTypeDonor == 3 ), pool.incompatiblePairs[i].donor_afam,\
        pool.incompatiblePairs[i].donor_age, pool.incompatiblePairs[i].donor_sex[0], pool.incompatiblePairs[i].donor_cig_use[0], \
        pool.incompatiblePairs[i].rec_sex[0], pool.incompatiblePairs[i].donor_weight, pool.incompatiblePairs[i].rec_weight, \
        pool.incompatiblePairs[i].donor_bmi)
    for j in range(K):
        if i == j:
            matches[i+T].append(0)
            continue
        compatible_1 = functions.are_blood_compatible(pool.incompatiblePairs[i].bloodTypeDonor, pool.incompatiblePairs[j].bloodTypePatient) \
                  and (not pool.incompatiblePairs[i].saidman.isPositiveCrossmatch(pool.incompatiblePairs[j].patientCPRA))
        compatible_2 = functions.are_blood_compatible(pool.incompatiblePairs[j].bloodTypeDonor, pool.incompatiblePairs[i].bloodTypePatient) \
                    and (not pool.incompatiblePairs[j].saidman.isPositiveCrossmatch(pool.incompatiblePairs[i].patientCPRA))
        if not (compatible_1 and compatible_2): 
            matches[i+T].append(0)
            continue
        donor_rec_abo_comp = functions.are_blood_compatible(pool.incompatiblePairs[i].bloodTypeDonor, pool.incompatiblePairs[j].bloodTypePatient)
        donor_rec_HLA_B_mis = gen.gen_donor_rec_HLA_B_mis(0)
        donor_rec_HLA_DR_mis = gen.gen_donor_rec_HLA_DR_mis(0)
        donor_rec_weight_ratio = util.calculate_weight_ratio(pool.incompatiblePairs[i].donor_weight, pool.incompatiblePairs[j].rec_weight)

        LKDPI_1 = util.calculate_lkdpi(pool.incompatiblePairs[i].donor_age, pool.incompatiblePairs[i].donor_afam, 
                                        pool.incompatiblePairs[i].donor_bmi, pool.incompatiblePairs[i].donor_cig_use,
                                          pool.incompatiblePairs[i].donor_sex, pool.incompatiblePairs[j].rec_sex,
                                          pool.incompatiblePairs[i].donor_sbp, donor_rec_abo_comp,
                                          0, pool.incompatiblePairs[i].donor_egfr, #assumed the donor and recip are NOT related
                                          donor_rec_HLA_B_mis, donor_rec_HLA_DR_mis,
                                          donor_rec_weight_ratio)
        donor_rec_abo_comp = functions.are_blood_compatible(pool.incompatiblePairs[j].bloodTypeDonor, pool.incompatiblePairs[i].bloodTypePatient)
        donor_rec_HLA_B_mis = gen.gen_donor_rec_HLA_B_mis(0)
        donor_rec_HLA_DR_mis = gen.gen_donor_rec_HLA_DR_mis(0)
        donor_rec_weight_ratio = util.calculate_weight_ratio(pool.incompatiblePairs[j].donor_weight, pool.incompatiblePairs[i].rec_weight)

        LKDPI_2 = util.calculate_lkdpi(pool.incompatiblePairs[j].donor_age, pool.incompatiblePairs[j].donor_afam, 
                                        pool.incompatiblePairs[j].donor_bmi, pool.incompatiblePairs[j].donor_cig_use,
                                          pool.incompatiblePairs[j].donor_sex, pool.incompatiblePairs[i].rec_sex,
                                          pool.incompatiblePairs[j].donor_sbp, donor_rec_abo_comp,
                                          0, pool.incompatiblePairs[j].donor_egfr, #assumed the donor and recip are NOT related
                                          donor_rec_HLA_B_mis, donor_rec_HLA_DR_mis,
                                          donor_rec_weight_ratio)
        matches[i+T].append(util.calculate_survival(LKDPI_1)+util.calculate_survival(LKDPI_2))


with open(filename, 'w') as f:
    f.write(json.dumps((K,T,matches, demo)))
