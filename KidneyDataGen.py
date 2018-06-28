#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Generates data for pools of incompatible pairs allowing up to 3 cycles. Outputs a json file with
the number of compatible pairs, number of incompatible pairs, 3D list of matches and demographic information 
about each pair.
"""
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
parser.add_argument('-o', '--output', default='data.dat')

args = parser.parse_args()

K = args.K
T = args.T
filename = args.output
pool = BJCSensitivityPool(T, K)
gen = DistributionGenerator()

misMatches = {}
positiveCrossMatches = {}
for i in range(T+K):
    for j in range(T+K):
        if j < T and i < T: continue
        misMatches[i+1,j+1] = (gen.gen_donor_rec_HLA_B_mis(0), gen.gen_donor_rec_HLA_DR_mis(0))


for t in pool.compatiblePairs:
    for i in pool.incompatiblePairs:
        positiveCrossMatches[t,i] = t.saidman.isPositiveCrossmatch(i.patientCPRA)
        positiveCrossMatches[i,t] = i.saidman.isPositiveCrossmatch(t.patientCPRA)
for i in pool.incompatiblePairs:
    for j in pool.incompatiblePairs:
        if i == j:
            positiveCrossMatches[i,j] = True
            continue
        positiveCrossMatches[i,j] = i.saidman.isPositiveCrossmatch(j.patientCPRA)

matchesUndirected = {}
matches2C = {}
matches3C = {}
demo = []

def generateDemo(pair):
    return (int(pair.bloodTypePatient == 0), int(pair.bloodTypePatient == 1 ), int(pair.bloodTypePatient == 2 ), \
         int(pair.bloodTypePatient == 3 ), int(pair.bloodTypeDonor == 0 ), int(pair.bloodTypeDonor == 1 ), \
         int(pair.bloodTypeDonor == 2 ), int(pair.bloodTypeDonor == 3 ), pair.donor_afam, pair.donor_age, \
         pair.donor_sex[0], pair.donor_cig_use[0], pair.rec_sex[0], pair.donor_weight, pair.rec_weight, \
         pair.donor_bmi, pair.donor_egfr, pair.donor_sbp, pair.patientCPRA)

def compatible(donor, recipient):
    return functions.are_blood_compatible(donor.bloodTypeDonor, recipient.bloodTypePatient) \
                  and (not positiveCrossMatches[donor, recipient])

#returns the EGS using LKDPI and assuming the donor and recipient are NOT related
def getLKDPI(donor, recipient, HLA_B_mis, HLA_DR_mis ):
    donor_rec_abo_comp = functions.are_blood_compatible(donor.bloodTypeDonor, recipient.bloodTypePatient)
    weight_ratio = util.calculate_weight_ratio(donor.donor_weight, recipient.rec_weight)
    return util.calculate_lkdpi(donor.donor_age, donor.donor_afam, donor.donor_bmi, donor.donor_cig_use, \
                                      donor.donor_sex, recipient.rec_sex, donor.donor_sbp, donor_rec_abo_comp, \
                                      0, donor.donor_egfr, HLA_B_mis, HLA_DR_mis, weight_ratio)
demo = [0]*(T+K)

#starting the cycle with a compatible pair
for i in range(T):
    matches2C[i+1,0] = util.calculate_survival(pool.compatiblePairs[i].LKDPI)
    matches3C[i+1,0,0] = util.calculate_survival(pool.compatiblePairs[i].LKDPI)
    matchesUndirected[i+1,0] = util.calculate_survival(pool.compatiblePairs[i].LKDPI)
    demo[i] = generateDemo(pool.compatiblePairs[i])
    for j in range(K):
        if compatible(pool.compatiblePairs[i], pool.incompatiblePairs[j]) and compatible(pool.incompatiblePairs[j], pool.compatiblePairs[i]): 
            lkdpi_ic = getLKDPI(pool.incompatiblePairs[j], pool.compatiblePairs[i], misMatches[j+T+1,i+1][0], misMatches[j+T+1,i+1][1])
            #compatibes will only donate with incentive (note: lower LKDPI is better)
            if lkdpi_ic < pool.compatiblePairs[i].LKDPI:
                lkdpi_ci = getLKDPI(pool.compatiblePairs[i], pool.incompatiblePairs[j], misMatches[i+1,j+T+1][0], misMatches[i+1,j+T+1][1])
                matches2C[i+1,j+1] = util.calculate_survival(lkdpi_ic) + util.calculate_survival(lkdpi_ci)
                matches3C[i+1,j+1,0] = util.calculate_survival(lkdpi_ic) + util.calculate_survival(lkdpi_ci)
                matchesUndirected[i+1,j+T+1] = util.calculate_survival(lkdpi_ci)
                matchesUndirected[j+T+1,i+1] = util.calculate_survival(lkdpi_ci)
        for k in range(K):
            if j == k: 
                continue
            #compatible[i] donates to incompatible[j] who donates to incompatible[k] who donates back to compatible[i]
            if compatible(pool.compatiblePairs[i], pool.incompatiblePairs[j]) and compatible(pool.incompatiblePairs[j], pool.incompatiblePairs[k]) \
            and compatible(pool.incompatiblePairs[k], pool.compatiblePairs[i]):
                lkdpi_ic = getLKDPI(pool.incompatiblePairs[k], pool.compatiblePairs[i], misMatches[k+T+1,i+1][0], misMatches[k+T+1,i+1][1])
                #if there is an incentive for the compatible pair to join the cycle
                if lkdpi_ic < pool.compatiblePairs[i].LKDPI:                    
                    lkdpi_ii = getLKDPI(pool.incompatiblePairs[j], pool.incompatiblePairs[k], misMatches[j+T+1,k+T+1][0], misMatches[j+T+1,k+T+1][1])
                    lkdpi_ci = getLKDPI(pool.compatiblePairs[i], pool.incompatiblePairs[j], misMatches[i+1,j+T+1][0], misMatches[i+1,j+T+1][1])
                    matches3C[i+1,j+1,k+1] = util.calculate_survival(lkdpi_ic) + util.calculate_survival(lkdpi_ii) + util.calculate_survival(lkdpi_ci)
            


for i in range(K):
    demo[i+T] = generateDemo(pool.incompatiblePairs[i])
    for j in range(K):
        if i == j: continue
        else:
            if compatible(pool.incompatiblePairs[i], pool.incompatiblePairs[j]) and compatible(pool.incompatiblePairs[j], pool.incompatiblePairs[i]):
                lkdpi_1 = getLKDPI(pool.incompatiblePairs[i], pool.incompatiblePairs[j], misMatches[i+T+1,j+T+1][0], misMatches[i+T+1,j+T+1][1])
                lkdpi_2 = getLKDPI(pool.incompatiblePairs[j], pool.incompatiblePairs[i], misMatches[j+T+1,i+T+1][0], misMatches[j+T+1,i+T+1][1])
                matches2C[i+T+1,j+1] = util.calculate_survival(lkdpi_1)+util.calculate_survival(lkdpi_2)
                matches3C[i+T+1,j+1,0] = util.calculate_survival(lkdpi_1)+util.calculate_survival(lkdpi_2)
                matchesUndirected[i+T+1,j+T+1] = util.calculate_survival(lkdpi_1)
                matchesUndirected[j+T+1,i+T+1] = util.calculate_survival(lkdpi_2)
            for k in range(K):
                if i == k or j == k: continue
                if compatible(pool.incompatiblePairs[i], pool.incompatiblePairs[j]) and compatible(pool.incompatiblePairs[j], pool.incompatiblePairs[k]) \
                and compatible(pool.incompatiblePairs[k], pool.compatiblePairs[i]):
                    lkdpi_1 = getLKDPI(pool.incompatiblePairs[i], pool.incompatiblePairs[j], misMatches[i+T+1,j+T+1][0], misMatches[i+T+1,j+T+1][1])
                    lkdpi_2 = getLKDPI(pool.incompatiblePairs[j], pool.incompatiblePairs[k], misMatches[j+T+1,k+T+1][0], misMatches[j+T+1,k+T+1][1])
                    lkdpi_3 = getLKDPI(pool.incompatiblePairs[k], pool.incompatiblePairs[i], misMatches[k+T+1,i+T+1][0], misMatches[k+T+1,i+T+1][1])
                    matches3C[i+T+1,j+1,k+1] = util.calculate_survival(lkdpi_1) + util.calculate_survival(lkdpi_2) + util.calculate_survival(lkdpi_3)

with open(filename, 'wb') as f:
    pickle.dump((K, T, matches2C, matches3C, demo, misMatches, matchesUndirected), f)
