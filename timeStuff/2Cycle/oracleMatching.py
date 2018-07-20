#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Takes a directory (default) or a single file of data and optimizes the model for count or quality.
"""
import json
import pickle
import argparse
import numpy as np
import pulp
from pulp import lpSum

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
    matches = d[3]
    directed_matches = d[7]
    T = num_compat
    model = pulp.LpProblem('oracle matching', pulp.LpMaximize)
    if args.incompatibleOnly:
        matchVars = [v for v in matches if v[0] > T]
    else:
        matchVars = [v for v in matches]
    x = pulp.LpVariable.dicts('match',matchVars,lowBound=0, upBound = 1, cat = pulp.LpInteger)
    matchVars = set(matchVars)
    if (args.quality):
        model += lpSum(x[v]*matches[v] for v in matchVars)
    else:
        model += lpSum(x[v]*COUNT(v) for v in matchVars)
    
    
    for t in range(1, num_pairs+1):
        model += lpSum(x[t,i] for i in range(num_incompat+1) if (t,i) in matchVars) <= 1,'only match with one '+str(t)
    
    for i in range(1,num_incompat+1):
        model += lpSum(x[t,i] for t in range(1,num_pairs+1) if (t,i) in matchVars) + \
                lpSum(x[i+T,j] for j in range(1,num_incompat+1) if (i+T,j) in matchVars) <= 1, 'symetry '+str(i)
  
    model.solve()
    
    used_incompat = set()      
    for v in matchVars:
        if round(x[v].value()) != 0:
            #if there is a compatible pair in the match
            if (v[0] <= T):
                #if compatible matched with itself
                if (v[1] == 0):
                    agentInfo += "C" + str(v[0]) + "\t" + str(0) + "\t" + str(directed_matches[v[0],0]) + "\t" \
                    + "C" + "\t" + str(directed_matches[v[0],0]) + "\t" + "C" + "\t" + str(demo[v[0]-1) + "\t" + str(0) + "\n"
                #compatible and incompatible
                else:
                    used_incompat.add(v[1])
                    agentInfo += "C" + str(v[0]) + "\t" + str(0) + "\t" + str(directed_matches[v[1]+T,v[0]]) + "\t" \
                    + "I" + "\t" + str(directed_matches[v[0],v[1]+T]) + "\t" + "I" + "\t" + str(demo[v[0]-1) + "\t" + str(0) + "\n"
                    agentInfo += "I" + str(v[1]) + "\t" + str(0) + "\t" + str(directed_matches[v[0],v[1]+T]) + "\t" \
                    + "C" + "\t" + str(directed_matches[v[1]+T,v[0]]) + "\t" + "C" + "\t" + str(demo[v[1]+C-1) + "\t" + str(departure_times[v[1]-1]) + "\n"
            #incompatible and incompatible
            else:
                used_incompat.add(v[0]-T)
                used_incompat.add(v[1])
                agentInfo += "I" + str(v[0]-T) + "\t" + str(0) + "\t" + str(directed_matches[v[1]+T,v[0]]) + "\t" \
                + "I" + "\t" + str(directed_matches[v[0],v[1]+T]) + "\t" + "I" + "\t" + str(demo[v[0]-1) + "\t" + str(departure_times[v[0]-1])) + "\n"
                agentInfo += "I" + str(v[1]) + "\t" + str(0) + "\t" + str(directed_matches[v[0],v[1]+T]) + "\t" \
                + "I" + "\t" + str(directed_matches[v[1]+T,v[0]]) + "\t" + "I" + "\t" + str(demo[v[1]+C-1) + "\t" + str(departure_times[v[1]-1]) + "\n"                
                
                
    if args.incompatibleOnly:
        quality = sum(x[v].value()*matches[v] for v in matchVars) + sum(matches[t,0] for t in range(1,T+1))
        count = sum(COUNT(v)*x[v].value() for v in matchVars) + T
    else:
        quality = sum(x[v].value()*matches[v] for v in matchVars)
        count = sum(COUNT(v)*x[v].value() for v in matchVars)
    
    for i in range(1,num_incompat+1):
        if i not in used_incompat:
            agentInfo += "I" + str(i) + "\t" + str(52) + "\t" + str(0) + "\t" \
                + "N" + "\t" + str(0) + "\t" + "N" + "\t" + str(demo[i+C-1) + "\t" + str(departure_times[i+C-1]) + "\n"
    
    results += str(count) + "\t" + str(quality) + "\n"

if args.output:
    with open(args.output, 'w') as f:
        f.write(results)
else:
    print results

if args.agents:
    with open(args.agents, 'w') as f:
        f.write(agentInfo)   
