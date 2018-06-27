#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Takes a directory (default) or a single file of data and optimizes the model for count or quality.
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
    with open (args.inputFile, 'rb') as f:
        data.append(pickle.load(f))
else:
    for fn in os.listdir(args.inputDir):
        if fn == '.DS_Store': continue
        with open(args.inputDir+'/'+fn, 'rb') as f:
            data.append(pickle.load(f))
def COUNT(v):
    if v[1] == 0:
        return 1
    if v[2] == 0:
        return 2
    return 3

pastData = []
for d in data:    
    num_incompat = d[0]
    num_compat = d[1]
    num_pairs = num_incompat + num_compat
    matches = d[2]
    T = num_compat
    model = Model('Kideny Optimizer')
    matchVars = {}
    for v in matches:
        matchVars[v] = model.addVar(vtype = GRB.CONTINUOUS, lb = 0, ub=1,  name = "match_" + str(v))
    
    model.addConstrs((quicksum(matchVars[t,i] for i in range(num_incompat+1) if (t,i) in matchVars) <= 1 \
                      for t in range(num_pairs)), "Only match with one pair")
    
    model.addConstrs((quicksum(matchVars[t,i] for t in range(num_pairs) if (t,i) in matchVars) + \
                      quicksum(matchVars[i+T-1,j] for j in range(1,num_incompat+1) if (i+T-1,j) in matchVars) <= 1 \
                      for i in range(1,num_incompat+1)), "undirected graph")
    if (args.quality):
        obj = quicksum(matchVars[v]*matches[v] for v in matchVars)
    else:
        obj = quicksum(COUNT(v)*matchVars[v] for v in matchVars)
  
    model.setObjective(obj, GRB.MAXIMIZE) 
    model.optimize()
    
    num_matches = 0
    for v in model.getVars():
        if v.X != 0:
            num_matches += 1
    quality = 0
    for i in range(num_pairs):
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