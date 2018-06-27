#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 27 09:07:27 2018

@author: kelseylieberman
"""

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Takes a list of .dat files from KidneyDataGen and optimizes the model for count or quality.
Returns count and quality for each population allowing cycles of up to three pairs.
"""
import pickle
import argparse
from gurobipy import *
import numpy as np

parser = argparse.ArgumentParser(description="Optimizes Kidney Exchange given by input file")
parser.add_argument('--inputFiles', nargs='+', default = ["data.dat"], help="list of .dat files to be used as input. List of number of \
                    incompatible pairs, number of compatible pairs, and list of all possible pairs with donor, recipient, egs")
parser.add_argument('--quality', action='store_true', help="Optimize for quality")
parser.add_argument('-o', '--output', default='data.dat')
args=parser.parse_args()

data = []
for fn in args.inputFiles:
    with open (fn, 'rb') as f:
        data.append(pickle.load(f))

def COUNT(v):
    if v[1] == 0:
        return 1
    if v[2] == 0:
        return 2
    return 3

results = ''
for d in data:    
    num_incompat = d[0]
    num_compat = d[1]
    num_pairs = num_incompat + num_compat
    matches = d[3]
    T = num_compat
    model = Model('Kideny Optimizer')
    matchVars = {}
    for v in matches:
        matchVars[v] = model.addVar(vtype = GRB.CONTINUOUS, lb = 0, ub=1,  name = "match_" + str(v))
    
    model.addConstrs((quicksum(matchVars[t,i,j] for i in range(num_incompat+1) for j in range(num_incompat+1) if (t,i,j) in matchVars) <= 1 \
                      for t in range(1,num_pairs+1)), "Only match with one pair")
    
    model.addConstrs((quicksum(matchVars[t,i,j] for t in range(1,num_pairs+1) for j in range(num_incompat+1) if (t,i,j) in matchVars) \
                     + quicksum(matchVars[t,j,i] for t in range(1,num_pairs+1) for j in range(1, num_incompat+1) if (t,j,i) in matchVars)\
                     + quicksum(matchVars[i+T,j,jp] for j in range(1,num_incompat+1) for jp in range(num_incompat+1) if (i+T,j,jp) in matchVars)\
                     <= 1 for i in range(1,num_incompat+1)), "undirected graph")

    if (args.quality):
        obj = quicksum(matchVars[v]*matches[v] for v in matchVars)

    else:
        obj = quicksum(COUNT(v)*matchVars[v] for v in matchVars)
        count = obj
        
    model.setObjective(obj, GRB.MAXIMIZE) 
    model.optimize()
    
    quality = sum(matchVars[v].X*matches[v] for v in matchVars)
    count = sum(COUNT(v)*matchVars[v].X for v in matchVars)
    
    results += str(quality) + "\t" + str(count) + "\n"

if args.output:
    with open(args.output, 'w') as f:
        f.write(results)
else:
    print results
