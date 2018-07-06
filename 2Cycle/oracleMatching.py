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
parser.add_argument('-o', '--output')
parser.add_argument('--incompatibleOnly', action = 'store_true', help = "Run the oracle on only the incompatible pairs")
parser.add_argument('--agents', help='output the quality of each agent to this file (.csv)')
args=parser.parse_args()

data = []
        
def COUNT(v):
    if v[1] == 0:
        return 1
    return 2

pastData = []
results = ''
agentInfo = ''
for fn in args.inputFiles:    
    with open(fn, 'rb') as f:
        d = pickle.load(f)
    num_incompat = d[0]
    num_compat = d[1]
    num_pairs = num_incompat + num_compat
    matches = d[2]
    directed_matches = d[6]
    T = num_compat
    model = Model('Kideny Optimizer')
    matchVars = {}
    for v in matches:
        if args.incompatibleOnly and v[0] <= T: continue
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
            
    for v in matchVars:
        if round(matchVars[v].X) != 0:
            #if there is a compatible pair in the match
            if (v[0] <= T):
                #if compatible matched with itself
                if (v[1] == 0):
                    agentInfo += "C" + str(v[0]) + "\t" + str(0) + "\t" + str(directed_matches[v[0],0]) + "\t" \
                    + "C" + "\t" + str(directed_matches[v[0],0]) + "\t" + "C" + "\n"
                #compatible and incompatible
                else:
                    agentInfo += "C" + str(v[0]) + "\t" + str(0) + "\t" + str(directed_matches[v[1]+T,v[0]]) + "\t" \
                    + "I" + "\t" + str(directed_matches[v[0],v[1]+T]) + "\t" + "I" + "\n"
                    agentInfo += "I" + str(v[1]) + "\t" + str(0) + "\t" + str(directed_matches[v[0],v[1]+T]) + "\t" \
                    + "C" + "\t" + str(directed_matches[v[1]+T,v[0]]) + "\t" + "C" + "\n"
            #incompatible and incompatible
            else:
                agentInfo += "I" + str(v[0]-T) + "\t" + str(0) + "\t" + str(directed_matches[v[1]+T,v[0]]) + "\t" \
                + "I" + "\t" + str(directed_matches[v[0],v[1]+T]) + "\t" + "I" + "\n"
                agentInfo += "I" + str(v[1]) + "\t" + str(0) + "\t" + str(directed_matches[v[0],v[1]+T]) + "\t" \
                + "I" + "\t" + str(directed_matches[v[1]+T,v[0]]) + "\t" + "I" + "\n"                
                
                
    if args.incompatibleOnly:
        quality = sum(matchVars[v].X*matches[v] for v in matchVars) + sum(matches[t,0] for t in range(1,T+1))
        count = sum(COUNT(v)*matchVars[v].X for v in matchVars) + T
    else:
        quality = sum(matchVars[v].X*matches[v] for v in matchVars)
        count = sum(COUNT(v)*matchVars[v].X for v in matchVars)
    
    results += str(count) + "\t" + str(quality) + "\n"

if args.output:
    with open(args.output, 'w') as f:
        f.write(results)
else:
    print results

if args.agents:
    with open(args.agents, 'w') as f:
        f.write(agentInfo)   
