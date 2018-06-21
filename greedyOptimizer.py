#!/usr/bin/env python2
# -*- coding: utf-8 -*-
import json
import argparse
from gurobipy import *
import numpy as np
from sets import Set

parser = argparse.ArgumentParser(description="Optimizes Kidney Exchange given by input file using a simple greedy algorithm")
parser.add_argument('--inputFile', nargs='?', help="JSON File to be used as input. List of number of \
                    incompatible pairs, number of compatible pairs, list of quality(egs) of all possible pairs \
                    and demographic information. File created in KidneyDataGen")
parser.add_argument('--quality', action='store_true', help="Optimize for quality")
parser.add_argument('-i', '--inputDir', nargs='?', default='data', help='input directory to look for data files')
parser.add_argument('-o', '--output', help='write results to this file')
args=parser.parse_args()

data = []
if args.inputFile:
    with open (args.inputFile, 'r') as f:
        data.append(json.load(f))
else:
    for fn in os.listdir(args.inputDir):
        if fn == '.DS_Store': continue
        with open(args.inputDir+'/'+fn, 'r') as f:
            data.append(json.load(f))

pastData = []

for d in data:    
    num_incompat = d[0]
    num_compat = d[1]
    num_pairs = num_incompat + num_compat
    matches = d[2]
    T = num_compat
    
    quality = 0
    num_compat_to_self = 0
    num_compat_to_incompat = 0
    num_incompat_to_compat = 0
    num_incompat_to_incompat = 0
    
    #greedy algorithm for compatible pairs by order of index
    used_incompat = Set([])
    for i in range(num_compat):
        max_index = 0
        for j in range(1,num_incompat+1):
            if matches[i][j] > matches[i][max_index] and j not in used_incompat:
                max_index = j
        quality += matches[i][max_index]
        used_incompat.add(max_index)
        if max_index == 0:
            num_compat_to_self += 1
        else:
            num_compat_to_incompat += 1
            num_incompat_to_compat += 1
    
    #gurobi optimization for remaining incompatible pairs
    model = Model('Kideny Optimizer')
    matchVars = {}
    incompat_matches = []
    for i in range(num_incompat):
        incompat_matches.append(matches[i + num_compat][1:])

    for i in range(num_incompat):
        if i in used_incompat: continue
        for j in range(num_incompat):
            if j not in used_incompat and incompat_matches[i][j] != 0:
                matchVars[(i,j)] = model.addVar(vtype = GRB.CONTINUOUS, lb = 0, ub=1,  name = "incompat_match_" + str((i,j)))
    
    model.addConstrs((quicksum(matchVars[i,j] for j in range(num_incompat) if (i,j) in matchVars) <= 1 \
                      for i in range(num_incompat)), "Only match with one pair")
    
    model.addConstrs((matchVars[(i,j)] == matchVars[(j,i)] for i in range(num_incompat) for j in range(num_incompat) \
                                if (i,j) in matchVars and (j,i) in matchVars), "undirected graph")

    if (args.quality):
        obj = quicksum(matchVars[i,j]*matches[i][j] for i in range(num_incompat) for j in range(num_incompat) if (i,j) in matchVars and (j,i) in matchVars)
    else:
        obj = quicksum(matchVars[i,j] for i in range(num_pairs) for j in range(num_incompat+1) if (i,j) in matchVars)
        
    model.setObjective(obj, GRB.MAXIMIZE) 
    model.optimize()
    
    quality += obj.getValue()
    for v in model.getVars():
        if v.X != 0:
            num_incompat_to_incompat += 1
            print v.varName
    num_matches = num_compat_to_self + num_compat_to_incompat + num_incompat_to_compat + num_incompat_to_incompat
    pastData.append((quality, num_matches))
    if args.output:
        with open(args.output, 'a') as f:
            f.write(str(quality) + "\t" + str(num_matches) + "\t" +  str(num_incompat_to_incompat) +"\n")
                     #num_compat_to_self, num_compat_to_incompat, num_incompat_to_compat, num_incompat_to_incompat))
avgs = np.mean(pastData, axis=0)
stdevs = np.std(pastData, axis=0)
s = ''
results = ''
for std in stdevs:
    s += str(std)+"\t"
s += "\n"
results = s+"\n\n\n"+results
s = ''
for a in avgs:
    s += str(a)+"\t"
s+="\n"
results = s+results

if args.output:
    with open(args.output, 'a') as f:
        f.write(results)
else:
    print results
