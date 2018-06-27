import argparse
from SaidmanCompatibleGenerator import SaidmanCompatibleGenerator
from SaidmanCompatibleGenerator import BJCSensitivityPool
from SaidmanCompatibleGenerator import BJCpool
from DistributionGenerator import *
import functions
import util
import pickle

parser = argparse.ArgumentParser(description="Generates Donor recipient pairs and a quality pool for optimization")
parser.add_argument('-K', '--num_incompatible', default=100, dest='K', type = int)
parser.add_argument('-T', '--num_compatible', default=100, dest='T', type = int)
parser.add_argument('-o', '--output', default='data.dat') #.dat because returns a binary file

args = parser.parse_args()

K = args.K
T = args.T
filename = args.output
pool = BJCSensitivityPool(T, K)
gen = DistributionGenerator()

def generateDemo(pair):
    return (int(pair.bloodTypePatient == 0), int(pair.bloodTypePatient == 1 ), int(pair.bloodTypePatient == 2 ), \
         int(pair.bloodTypePatient == 3 ), int(pair.bloodTypeDonor == 0 ), int(pair.bloodTypeDonor == 1 ), \
         int(pair.bloodTypeDonor == 2 ), int(pair.bloodTypeDonor == 3 ), pair.donor_afam, pair.donor_age, \
         pair.donor_sex[0], pair.donor_cig_use[0], pair.rec_sex[0], pair.donor_weight, pair.rec_weight, \
         pair.donor_bmi, pair.donor_egfr, pair.donor_sbp, pair.patientCPRA)

def compatible(donor, recipient):
    return functions.are_blood_compatible(donor.bloodTypeDonor, recipient.bloodTypePatient) \
                  and (not donor.saidman.isPositiveCrossmatch(recipient.patientCPRA))

#returns the EGS using LKDPI and assuming the donor and recipient are NOT related
def getLKDPI(donor, recipient):
    donor_rec_abo_comp = functions.are_blood_compatible(donor.bloodTypeDonor, recipient.bloodTypePatient)
    donor_rec_HLA_B_mis = gen.gen_donor_rec_HLA_B_mis(0)
    donor_rec_HLA_DR_mis = gen.gen_donor_rec_HLA_DR_mis(0)
    weight_ratio = util.calculate_weight_ratio(donor.donor_weight, recipient.rec_weight)
    return util.calculate_lkdpi(donor.donor_age, donor.donor_afam, donor.donor_bmi, donor.donor_cig_use, \
                                      donor.donor_sex, recipient.rec_sex, donor.donor_sbp, donor_rec_abo_comp, \
                                      0, donor.donor_egfr, donor_rec_HLA_B_mis, donor_rec_HLA_DR_mis, weight_ratio)
matches = {}
demo = [0]*(T+K)

for i in range(T):
    matches[i,0] = util.calculate_survival(pool.compatiblePairs[i].LKDPI)
    demo[i] = generateDemo(pool.compatiblePairs[i])
    for j in range(K):
        if compatible(pool.compatiblePairs[i], pool.incompatiblePairs[j]) and compatible(pool.incompatiblePairs[j], pool.compatiblePairs[i]): 
            lkdpi_ic = getLKDPI(pool.incompatiblePairs[j], pool.compatiblePairs[i])
            #compatibes will only donate with incentive (note: lower LKDPI is better)
            if lkdpi_ic < pool.compatiblePairs[i].LKDPI:
                lkdpi_ci = getLKDPI(pool.compatiblePairs[i], pool.incompatiblePairs[j])
                matches[i,j+1]= util.calculate_survival(lkdpi_ic) + util.calculate_survival(lkdpi_ci)
            

for i in range(K):
    matches[i+T][0] = 0
    demo[i+T] = generateDemo(pool.incompatiblePairs[i])
    for j in range(K):
        if i == j:
            continue
        else:
            if compatible(pool.incompatiblePairs[i], pool.incompatiblePairs[j]) and compatible(pool.incompatiblePairs[j], pool.incompatiblePairs[i]):
                lkdpi_1 = getLKDPI(pool.incompatiblePairs[i], pool.incompatiblePairs[j])
                lkdpi_2 = getLKDPI(pool.incompatiblePairs[j], pool.incompatiblePairs[i])
                matches[i+T,j+1] = util.calculate_survival(lkdpi_1)+util.calculate_survival(lkdpi_2)

with open(filename, 'wb') as f:
    pickle.dump((K, T, matches, demo), f)