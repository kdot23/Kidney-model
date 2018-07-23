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
from SaidmanCompatibleGenerator import BJCPair
from DistributionGenerator import *
import functions
import util
import pickle
import random
import numpy as np

parser = argparse.ArgumentParser(description="Generates Donor recipient pairs and a quality pool for optimization")
parser.add_argument('-K', '--num_incompatible', default=100, dest='K', type = int)
parser.add_argument('-T', '--time', default=100, dest='T', type = int)
parser.add_argument('-o', '--output', default='data.dat')
parser.add_argument('--l1', default=1., type=float, help = "Incompatible arrival rate")
parser.add_argument('--l2', default=1., type=float, help = "Compatible arrival rate")
parser.add_argument('-L', '--mean_life', type=int, default=75, help='Mean lifespan of an incompatible pair')
parser.add_argument('-S', '--stdev_life', type=int, default=5, help='Stdev of lifespans of incompatible pairs')

args = parser.parse_args()

K = args.K
T = args.T
filename = args.output
pool = BJCSensitivityPool(0, K)
gen = DistributionGenerator()

incompatibles = pool.incompatiblePairs
compatibles = []
departureTimesIncompat = []

for t in range(1,T+1):
    incompatArrivals = np.random.poisson(args.l1)
    for _ in range(incompatArrivals):
        pair = BJCPair(arrival_time=t)
        while pair.compatible:
            pair = BJCPair(arrival_time=t)
        incompatibles.append(pair)
    compatArrivals = np.random.poisson(args.l2)
    for _ in range(compatArrivals):
        pair = BJCPair(arrival_time=t)
        while not pair.compatible:
            pair = BJCPair(arrival_time=t)
        compatibles.append(pair)

for i in range(len(incompatibles)):
    time = round(max(0,np.random.normal(loc=args.mean_life, scale=args.stdev_life)))
    incompatibles[i].life_span = time
    departureTimesIncompat.append(incompatibles[i].arrival_time+time)


C = len(compatibles)
I = len(incompatibles)
misMatches = {}
positiveCrossMatches = {}
for i in range(len(compatibles)+len(incompatibles)):
    for j in range(len(compatibles)+len(incompatibles)):
        if j < C and i < C: continue
        misMatches[i+1,j+1] = (gen.gen_donor_rec_HLA_B_mis(0), gen.gen_donor_rec_HLA_DR_mis(0))


for t in compatibles:
    for i in incompatibles:
        positiveCrossMatches[t,i] = t.saidman.isPositiveCrossmatch(i.pr_PraIncompatiblity)
        positiveCrossMatches[i,t] = i.saidman.isPositiveCrossmatch(t.pr_PraIncompatiblity)
for i in incompatibles:
    for j in incompatibles:
        if i == j:
            positiveCrossMatches[i,j] = True
            continue
        positiveCrossMatches[i,j] = i.saidman.isPositiveCrossmatch(j.pr_PraIncompatiblity)

matchesDirected = {}
matches2C = {}
matches3C = {}
demo = []

def generateDemo(pair):
    return (int(pair.bloodTypePatient == 0), int(pair.bloodTypePatient == 1 ), int(pair.bloodTypePatient == 2 ), \
         int(pair.bloodTypePatient == 3 ), int(pair.bloodTypeDonor == 0 ), int(pair.bloodTypeDonor == 1 ), \
         int(pair.bloodTypeDonor == 2 ), int(pair.bloodTypeDonor == 3 ), pair.donor_afam, pair.donor_age, \
         pair.donor_sex[0], pair.donor_cig_use[0], pair.rec_sex[0], pair.donor_weight, pair.rec_weight, \
         pair.donor_bmi, pair.donor_egfr, pair.donor_sbp, pair.pr_PraIncompatiblity, pair.life_span, pair.arrival_time)

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
demo = [0]*(C+I)

#starting the cycle with a compatible pair
for i in range(len(compatibles)):
    matches2C[i+1,0] = util.calculate_survival(compatibles[i].LKDPI)
    matches3C[i+1,0,0] = util.calculate_survival(compatibles[i].LKDPI)
    matchesDirected[i+1,0] = util.calculate_survival(compatibles[i].LKDPI)
    demo[i] = generateDemo(compatibles[i])
    for j in range(len(incompatibles)):
        if departureTimesIncompat[j] < compatibles[i].arrival_time or incompatibles[j].arrival_time > compatibles[i].arrival_time:
            continue
        if compatible(compatibles[i], incompatibles[j]):
            matchesDirected[i+1,j+C+1] = util.calculate_survival(getLKDPI(compatibles[i], incompatibles[j], misMatches[i+1,j+C+1][0], misMatches[i+1,j+C+1][1]))
        if compatible(incompatibles[j], compatibles[i]):
            matchesDirected[j+C+1,i+1] = util.calculate_survival(getLKDPI(incompatibles[j], compatibles[i], misMatches[j+C+1,i+1][0], misMatches[j+C+1,i+1][1]))
        if compatible(compatibles[i], incompatibles[j]) and compatible(incompatibles[j],compatibles[i]): 
            lkdpi_ic = getLKDPI(incompatibles[j], compatibles[i], misMatches[j+C+1,i+1][0], misMatches[j+C+1,i+1][1])
            #compatibes will only donate with incentive (note: lower LKDPI is better)
            if lkdpi_ic < compatibles[i].LKDPI:
                lkdpi_ci = getLKDPI(compatibles[i], incompatibles[j], misMatches[i+1,j+C+1][0], misMatches[i+1,j+C+1][1])
                matches2C[i+1,j+1] = util.calculate_survival(lkdpi_ic) + util.calculate_survival(lkdpi_ci)
                matches3C[i+1,j+1,0] = util.calculate_survival(lkdpi_ic) + util.calculate_survival(lkdpi_ci)
        for k in range(len(incompatibles)):
            if departureTimesIncompat[k] < compatibles[i].arrival_time or incompatibles[k].arrival_time > compatibles[i].arrival_time:
                continue
            if departureTimesIncompat[k] < incompatibles[j].arrival_time or incompatibles[k].arrival_time > departureTimesIncompat[j]:
                continue
            if j == k: 
                continue
            #compatible[i] donates to incompatible[j] who donates to incompatible[k] who donates back to compatible[i]
            if compatible(compatibles[i], incompatibles[j]) and compatible(incompatibles[j], incompatibles[k]) \
            and compatible(incompatibles[k], compatibles[i]):
                lkdpi_ic = getLKDPI(incompatibles[k], compatibles[i], misMatches[k+C+1,i+1][0], misMatches[k+C+1,i+1][1])
                #if there is an incentive for the compatible pair to join the cycle
                if lkdpi_ic < compatibles[i].LKDPI:                    
                    lkdpi_ii = getLKDPI(incompatibles[j], incompatibles[k], misMatches[j+C+1,k+C+1][0], misMatches[j+C+1,k+C+1][1])
                    lkdpi_ci = getLKDPI(compatibles[i], incompatibles[j], misMatches[i+1,j+C+1][0], misMatches[i+1,j+C+1][1])
                    matches3C[i+1,j+1,k+1] = util.calculate_survival(lkdpi_ic) + util.calculate_survival(lkdpi_ii) + util.calculate_survival(lkdpi_ci)
            


for i in range(len(incompatibles)):
    demo[i+C] = generateDemo(pool.incompatiblePairs[i])
    for j in range(len(incompatibles)):
        if departureTimesIncompat[j] < incompatibles[i].arrival_time or incompatibles[j].arrival_time > departureTimesIncompat[i]:
            continue
        if i == j: continue
        if compatible(incompatibles[i], incompatibles[j]):
            matchesDirected[i+C+1,j+C+1] = util.calculate_survival(getLKDPI(incompatibles[i],incompatibles[j], misMatches[i+C+1,j+C+1][0], misMatches[i+C+1,j+C+1][1]))
        if compatible(incompatibles[j], incompatibles[i]):
            matchesDirected[j+C+1,i+C+1] = util.calculate_survival(getLKDPI(incompatibles[j],incompatibles[i], misMatches[j+C+1,i+C+1][0], misMatches[j+C+1,i+C+1][1]))
        if compatible(incompatibles[i], incompatibles[j]) and compatible(incompatibles[j], incompatibles[i]):
            lkdpi_1 = getLKDPI(incompatibles[i], incompatibles[j], misMatches[i+C+1,j+C+1][0], misMatches[i+C+1,j+C+1][1])
            lkdpi_2 = getLKDPI(incompatibles[j], incompatibles[i], misMatches[j+C+1,i+C+1][0], misMatches[j+C+1,i+C+1][1])
            matches2C[i+C+1,j+1] = util.calculate_survival(lkdpi_1)+util.calculate_survival(lkdpi_2)
            matches3C[i+C+1,j+1,0] = util.calculate_survival(lkdpi_1)+util.calculate_survival(lkdpi_2)
        for k in range(len(incompatibles)):
            if departureTimesIncompat[k] < incompatibles[i].arrival_time or incompatibles[k].arrival_time > departureTimesIncompat[i]:
                continue
            if departureTimesIncompat[k] < incompatibles[j].arrival_time or incompatibles[k].arrival_time > departureTimesIncompat[j]:
                continue
            if i == k or j == k: continue
            if compatible(incompatibles[i], incompatibles[j]) and compatible(incompatibles[j], incompatibles[k]) \
            and compatible(incompatibles[k], incompatibles[i]):
                lkdpi_1 = getLKDPI(incompatibles[i], incompatibles[j], misMatches[i+C+1,j+C+1][0], misMatches[i+C+1,j+C+1][1])
                lkdpi_2 = getLKDPI(incompatibles[j], incompatibles[k], misMatches[j+C+1,k+C+1][0], misMatches[j+C+1,k+C+1][1])
                lkdpi_3 = getLKDPI(incompatibles[k], incompatibles[i], misMatches[k+C+1,i+C+1][0], misMatches[k+C+1,i+C+1][1])
                matches3C[i+C+1,j+1,k+1] = util.calculate_survival(lkdpi_1) + util.calculate_survival(lkdpi_2) + util.calculate_survival(lkdpi_3)

with open(filename, 'wb') as f:
    pickle.dump((I, C, T, matches2C, matches3C, demo, misMatches, matchesDirected,departureTimesIncompat), f)
