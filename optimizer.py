#import modules
import numpy as np
import os
import json
from gurobipy import *
import argparse
from SaidmanCompatibleGenerator import SaidmanCompatibleGenerator
from SaidmanCompatibleGenerator import BJCSensitivityPool
from SaidmanCompatibleGenerator import BJCpool
from DistributionGenerator import *
import functions
import util
import pickle

parser = argparse.ArgumentParser(description="Optimizes Kidney Exchange given by input file")
"""
parser.add_argument('--inputFile', nargs='?', help="JSON File to be used")
parser.add_argument('-i', '--inputDir', nargs='?', default='data', help='input directory to look for data files')
parser.add_argument('-o', '--outputFile', nargs='?', help="CSV File for results to be saved to")
parser.add_argument('--quality', action='store_true', help="Optimize for quality")

#BJCSensitivityPool 
pool = BJCSensitivityPool(100,100)
for i in pool.compatibleIDs:
   print pool.compatibleIDs[i]
"""
pool = BJCSensitivityPool(50,50)
print len(pool.compatibleIDs)
print len(pool.incompatibleIDs)
for i in range(len(pool.compatibleIDs)):
   print pool.compatibleIDs[i]
    
"""
args=parser.parse_args()
results = "#nm\tnum_C\tnum_CI\tnum_II\tq\tc_q\tci_q\tii_q\n"
if args.inputFile:
    with open(args.inputFile,'r') as f:
        data = [json.load(f)]
else:
    data = []
    for fn in os.listdir(args.inputDir):
        if fn == '.DS_Store': continue
        with open(args.inputDir+'/'+fn, 'r') as f:
            data.append(json.load(f))

pastData = []
for d in data:
    compat = d[0]
    incompat = d[1]
    Igraph = d[2]
    CIgraph = d[3]
    T = len(compat)
""" 
model = Model('Kideny Optimizer')
iMatches = {}
ciMatches = {}
iiLKDPI = {}
ciLKDPI = {}
icLKDPI = {}
cLKDPI = {}
#only create variables for edges which exist


for i in range(len(pool.incompatiblePairs)):
    for j in range(len(pool.incompatiblePairs)):
        if j > i:
            #possible += 1
            compatible_1 = functions.are_blood_compatible(pool.incompatiblePairs[i].bloodTypeDonor, pool.incompatiblePairs[j].bloodTypePatient) \
                      and (not pool.incompatiblePairs[i].saidman.isPositiveCrossmatch(pool.incompatiblePairs[j].patientCPRA))
            compatible_2 = functions.are_blood_compatible(pool.incompatiblePairs[j].bloodTypeDonor, pool.incompatiblePairs[i].bloodTypePatient) \
                        and (not pool.incompatiblePairs[j].saidman.isPositiveCrossmatch(pool.incompatiblePairs[i].patientCPRA))
            if compatible_1 and compatible_2:
                iMatches[(i,j)] = model.addVar(vtype = GRB.CONTINUOUS, lb = 0, ub = 1, name = "I" + str((i,j)))
            
            donor_rec_HLA_B_mis = gen.gen_donor_rec_HLA_B_mis(0)
            donor_rec_HLA_DR_mis = gen.gen_donor_rec_HLA_DR_mis(0)
            donor_rec_abo_comp_ij = functions.are_blood_compatible(pool.incompatiblePairs[i].bloodTypeDonor, pool.compatiblePairs[j].bloodTypePatient)
            donor_rec_weight_ratio_ij = util.calculate_weight_ratio(pool.incompatiblePairs[i].donor_weight, pool.incompatiblePairs[j].rec_weight)   
            LKDPI_ij = util.calculate_lkdpi(pool.incompatiblePairs[i].donor_age, pool.incompatiblePairs[i].donor_afam, 
                                            pool.incompatiblePairs[i].donor_bmi, pool.incompatiblePairs[i].donor_cig_use,
                                              pool.incompatiblePairs[i].donor_sex, pool.incompatiblePairs[j].rec_sex,
                                              pool.incompatiblePairs[i].donor_sbp, donor_rec_abo_comp_ij,
                                              0, pool.incompatiblePairs[i].donor_egfr, #assumed the donor and recip are NOT related
                                              donor_rec_HLA_B_mis, donor_rec_HLA_DR_mis,
                                              donor_rec_weight_ratio_ij)
            donor_rec_abo_comp_ji = functions.are_blood_compatible(pool.incompatiblePairs[j].bloodTypeDonor, pool.compatiblePairs[i].bloodTypePatient)
            donor_rec_weight_ratio_ji = util.calculate_weight_ratio(pool.incompatiblePairs[j].donor_weight, pool.incompatiblePairs[i].rec_weight)   
            LKDPI_ji = util.calculate_lkdpi(pool.incompatiblePairs[j].donor_age, pool.incompatiblePairs[j].donor_afam, 
                                            pool.incompatiblePairs[j].donor_bmi, pool.incompatiblePairs[j].donor_cig_use,
                                              pool.incompatiblePairs[j].donor_sex, pool.incompatiblePairs[i].rec_sex,
                                              pool.incompatiblePairs[j].donor_sbp, donor_rec_abo_comp_ji,
                                              0, pool.incompatiblePairs[j].donor_egfr, #assumed the donor and recip are NOT related
                                              donor_rec_HLA_B_mis, donor_rec_HLA_DR_mis,
                                              donor_rec_weight_ratio_ji)
            iiLKDPI[(i,j)] = LKDPI_ij
            iiLKDPI[(j,i)] = LKDPI_ji
                

gen = DistributionGenerator()

for i in range(len(pool.compatiblePairs)):
    for j in range(len(pool.incompatiblePairs)):
        compatible_1 = functions.are_blood_compatible(pool.compatiblePairs[i].bloodTypeDonor, pool.incompatiblePairs[j].bloodTypePatient) \
                  and (not pool.compatiblePairs[i].saidman.isPositiveCrossmatch(pool.incompatiblePairs[j].patientCPRA))
        compatible_2 = functions.are_blood_compatible(pool.incompatiblePairs[j].bloodTypeDonor, pool.compatiblePairs[i].bloodTypePatient) \
                    and (not pool.incompatiblePairs[j].saidman.isPositiveCrossmatch(pool.compatiblePairs[i].patientCPRA))
        if not (compatible_1 and compatible_2): continue
        #LKPDI of incompatible < LKPDI of compatible for variable to exist
        #possible += 1
        
        donor_rec_HLA_B_mis = gen.gen_donor_rec_HLA_B_mis(0)
        donor_rec_HLA_DR_mis = gen.gen_donor_rec_HLA_DR_mis(0)
        donor_rec_abo_comp_ci = functions.are_blood_compatible(pool.incompatiblePairs[j].bloodTypeDonor, pool.compatiblePairs[i].bloodTypePatient)
        donor_rec_weight_ratio_ci = util.calculate_weight_ratio(pool.incompatiblePairs[j].donor_weight, pool.compatiblePairs[i].rec_weight)
        LKDPI_CI = util.calculate_lkdpi(pool.incompatiblePairs[j].donor_age, pool.incompatiblePairs[j].donor_afam, 
                                        pool.incompatiblePairs[j].donor_bmi, pool.incompatiblePairs[j].donor_cig_use,
                                          pool.incompatiblePairs[j].donor_sex, pool.compatiblePairs[i].rec_sex,
                                          pool.incompatiblePairs[j].donor_sbp, donor_rec_abo_comp_ci,
                                          0, pool.incompatiblePairs[j].donor_egfr, #assumed the donor and recip are NOT related
                                          donor_rec_HLA_B_mis, donor_rec_HLA_DR_mis,
                                          donor_rec_weight_ratio_ci)
        
        if LKDPI_CI < pool.compatiblePairs[i].LKDPI:
            #compatible_pairs += 1
            ciMatches[(i,j)] = model.addVar(vtype = GRB.CONTINUOUS, lb = 0, ub = 1, name = "CI" + str((i,j)))
            donor_rec_abo_comp_ic = functions.are_blood_compatible(pool.compatiblePairs[i].bloodTypeDonor, pool.incompatiblePairs[j].bloodTypePatient)
            donor_rec_weight_ratio_ic = util.calculate_weight_ratio(pool.compatiblePairs[i].donor_weight, pool.incompatiblePairs[j].rec_weight)
            LKDPI_IC = util.calculate_lkdpi(pool.compatiblePairs[i].donor_age, pool.compatiblePairs[i].donor_afam, 
                                            pool.compatiblePairs[i].donor_bmi, pool.compatiblePairs[i].donor_cig_use,
                                              pool.compatiblePairs[i].donor_sex, pool.incompatiblePairs[j].rec_sex,
                                              pool.compatiblePairs[i].donor_sbp, donor_rec_abo_comp_ic,
                                              0, pool.compatiblePairs[i].donor_egfr, #assumed the donor and recip are NOT related
                                              donor_rec_HLA_B_mis, donor_rec_HLA_DR_mis,
                                              donor_rec_weight_ratio_ic)
            ciLKDPI[(j,i)] = LKDPI_CI           
            icLKDPI[(i,j)] = LKDPI_IC
            
    ciMatches[(i,-1)] = model.addVar(vtype = GRB.CONTINUOUS, lb = 0, ub = 1, name = "CI" + str((i,-1)))
    cLKDPI[i] = pool.compatiblePairs[i].LKDPI
    
    #compatible pair can only match with up to one incompatible pair or themselves
    model.addConstrs((quicksum(ciMatches[i,j] for j in range(-1,len(pool.incompatiblePairs)) if (i,j) in ciMatches)  <= 1 for i in range(len(pool.compatiblePairs))), "Compat  Matches")
    
    #incompatible pair can only match with up to one incompatible pair or one compatible pair
    model.addConstrs((quicksum(iMatches[i,j] for j in range(len(pool.incompatiblePairs)) if (i,j) in iMatches) + quicksum(ciMatches[j,i] for j in range(len(pool.incompatiblePairs)) if (j,i) in ciMatches) <= 1 for i in range(len(pool.compatiblePairs))), "Incompat Matches")
    
    # model.addConstrs((iMatches[i,j] == iMatches[j,i] for i in range(T) for j in range(T) if (i,j) in iMatches), "Symmetry") 
    #obj = quicksum(ciMatches[i,j] for k in ciMatches) + quicksum(iMatches[i,j] for k in iMatches) 
    obj = quicksum(ciMatches[k]*(icLKDPI[(k)] + ciLKDPI[(k[1],k[0])]) for k in ciMatches) 
    #+ quicksum(iMatches[k]*iiLKDPI[k] for k in iMatches) + quicksum(ciMatches([i,-1])*cLKDPI[i])
    model.setObjective(obj, GRB.MAXIMIZE) 
    model.optimize()
    print obj.getValue()
"""   
    if (args.quality):
        obj = quicksum((compat[i][2]+incompat[j][2])*ciMatches[i,j] for i in range(T) for j in range(T) if (i,j) in ciMatches) + quicksum(compat[i][2]*ciMatches[i,-1] for i in range(T) ) + quicksum(incompat[i][2]*iMatches[i,j] for i in range(T) for j in range(T) if (i,j) in iMatches)
    else:
        obj = quicksum(ciMatches[i,j] for i in range(len(pool.compatiblePairs)) for j in range(pool.incompatiblePairs) if (i,j) in ciMatches) + quicksum(ciMatches[i,-1] for i in range(len(pool.compatiblePairs)) ) + quicksum(iMatches[i,j] for i in range(len(incompatiblePairs)) for j in range(len(incompatiblePairs)) if (i,j) in iMatches)
    model.setObjective(obj, GRB.MAXIMIZE) 
    model.optimize()
    num_matches = 0
    quality = 0.
    c_quality = 0.
    ci_quality = 0.
    ii_quality = 0.
    num_II = 0
    num_CI = 0
    num_C = 0
    for v in model.getVars():
        if v.X != 0:
            #print str(v.Varname) + ' ' + str(v.X)
            num_matches += 1
    
    for var in iMatches:
        if (iMatches[var].X != 0):
            ii_quality += incompat[var[0]][2]
            num_II += 1
    for var in ciMatches:
        if (ciMatches[var].X != 0):
            if var[1] == -1:
                num_C += 1
                c_quality += compat[var[0]][2]
            else:
                num_CI += 1
                ci_quality += compat[var[0]][2] + incompat[var[1]][2]
    
    print num_matches
    print num_C
    print num_CI
    print num_II
    
    quality = c_quality + ci_quality + ii_quality
    print quality
    if num_C > 0:
        print c_quality/num_C
    if num_CI > 0:
        print ci_quality/num_CI
    if (num_II > 0):
        print ii_quality/num_II
    results += str(num_matches) + "\t" + str(num_C) + "\t" + str(num_CI) + "\t" + str(num_II) + "\t" + str(quality) + "\t" 
    if num_C > 0:
        results += str(c_quality/num_C) + "\t" 
    else:
        results += "na\t"
    if num_CI > 0:
        results += str(ci_quality/num_CI) + "\t"
    else:
        results += "na\t"
    if (num_II > 0):
        results += str(ii_quality/num_II)
    else:
        results += "na"
    results += "\n"
    pastData.append((num_matches,num_C, num_CI, num_II, quality, 
        c_quality/num_C if num_C else 0, ci_quality/num_CI if num_CI else 0, ii_quality/num_II if num_II else 0))

avgs = np.mean(pastData, axis=0)
stdevs = np.std(pastData, axis=0)

s = ''
for std in stdevs:
    s += str(std)+"\t"
s += "\n"
results = s+"\n\n\n"+results
s = ''
for a in avgs:
    s += str(a)+"\t"
s+="\n"
results = s+results
    
if args.outputFile:
    with open(args.outputFile,'w') as f:
        f.write(results)
else:
    print results
"""
