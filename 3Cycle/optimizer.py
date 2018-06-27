#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 27 09:07:27 2018

@author: kelseylieberman
"""

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Takes a directory (default) or a single file of data and optimizes the model for count or quality.
Allows cycles of up to three pairs.
"""
import json
import argparse
from gurobipy import *
import numpy as np

parser = argparse.ArgumentParser(description="Optimizes Kidney Exchange given by input file")
parser.add_argument('--inputFile', nargs='?', help="JSON File to be used as input. List of number of \
                    incompatible pairs, number of compatible pairs, and list of all possible pairs with donor, recipient, egs")
parser.add_argument('--quality', action='store_true', help="Optimize for quality")
parser.add_argument('-i', '--inputDir', nargs='?', default='data', help='input directory to look for data files')
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
    model = Model('Kideny Optimizer')
    matchVars = {}
    for t in range(num_pairs):
        for i in range(num_incompat+1):
            for j in range(num_incompat + 1):                
                if matches[t][i][j] != 0:
                    matchVars[(t,i,j)] = model.addVar(vtype = GRB.CONTINUOUS, lb = 0, ub=1,  name = "match_" + str((t,i,j)))
    
    model.addConstrs((quicksum(matchVars[t,i,j] for i in range(num_incompat+1) for j in range(num_incompat+1) if (t,i,j) in matchVars) <= 1 \
                      for t in range(num_pairs)), "Only match with one pair")
    
    model.addConstrs((quicksum(matchVars[t,i,j] for t in range(num_pairs) for j in range(num_incompat+1) if (t,i,j) in matchVars) \
                     + quicksum(matchVars[t,j,i] for t in range(num_pairs) for j in range(1, num_incompat+1) if (t,j,i) in matchVars)\
                     + quicksum(matchVars[i+T-1,j,jp] for j in range(1,num_incompat+1) for jp in range(num_incompat+1) if (i+T-1,j,jp) in matchVars)\
                     <= 1 for i in range(1,num_incompat+1)), "undirected graph")

    if (args.quality):
        obj = quicksum(matchVars[t,i,j]*matches[t][i][t] for t in range(num_pairs) for i in range(num_incompat+1) for j in range(num_incompat+1)\
                       if (t,i,j) in matchVars)
    else:
        obj = quicksum(matchVars[t,i,j] for t in range(num_pairs) for i in range(num_incompat+1) for j in range(num_incompat+1) if (t,i,j) in matchVars)
        
    model.setObjective(obj, GRB.MAXIMIZE) 
    model.optimize()
    
    num_matches = 0
    for v in model.getVars():
        if v.X != 0:
            num_matches += 1
    quality = 0
    for t in range(num_pairs):
        for i in range(num_incompat+1):
            for j in range(num_incompat + 1):
                if (i,j) in matchVars:
                    quality += matchVars[i,j].X*matches[i][j]
    pastData.append((obj.getValue(), quality, num_matches, quality/num_matches))
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
print results