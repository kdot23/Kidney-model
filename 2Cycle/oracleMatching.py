#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Takes a directory (default) or a single file of data and optimizes the model for count or quality.
"""
import json
import pickle
import argparse
from gurobipy import *
import numpy as np

parser = argparse.ArgumentParser(description="Optimizes Kidney Exchange given by input file")
parser.add_argument('--inputFiles', nargs='+', default = ["data.dat"], help="list of .dat files to be used as input. List of number of \
                    incompatible pairs, number of compatible pairs, and list of all possible pairs with donor, recipient, egs")
parser.add_argument('--quality', action='store_true', help="Optimize for quality")
parser.add_argument('-o', '--output', default='data.csv')
args=parser.parse_args()

data = []
        
def COUNT(v):
    if v[1] == 0:
        return 1
    return 2

pastData = []
results = ''
for fn in args.inputFiles:    
    with open(fn, 'rb') as f:
        d = pickle.load(f)
    num_incompat = d[0]
    num_compat = d[1]
    num_pairs = num_incompat + num_compat
    matches = d[2]
    T = num_compat
    model = Model('Kideny Optimizer')
    matchVars = {}
    for v in matches:
        matchVars[v] = model.addVar(vtype = GRB.BINARY, lb = 0, ub=1,  name = "match_" + str(v))
    
    model.addConstrs((quicksum(matchVars[t,i] for i in range(num_incompat+1) if (t,i) in matchVars) <= 1 \
                      for t in range(1,num_pairs+1)), "Only match with one pair")
    
    model.addConstrs((quicksum(matchVars[t,i] for t in range(1,num_pairs+1) if (t,i) in matchVars) + \
                      quicksum(matchVars[i+T,j] for j in range(1,num_incompat+1) if (i+T,j) in matchVars) <= 1 \
                      for i in range(1,num_incompat+1)), "undirected graph")
    if (args.quality):
        obj = quicksum(matchVars[v]*matches[v] for v in matchVars)
    else:
        obj = quicksum(COUNT(v)*matchVars[v] for v in matchVars)
  
    model.setObjective(obj, GRB.MAXIMIZE) 
    model.optimize()
            
    quality = sum(matchVars[v].X*matches[v] for v in matchVars)
    count = sum(COUNT(v)*matchVars[v].X for v in matchVars)
    
    results += str(count) + "\t" + str(quality) + "\n"

if args.output:
    with open(args.output, 'w') as f:
        f.write(results)
else:
    print results
