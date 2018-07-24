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
import numpy as np
import pulp
from pulp import lpSum
from zipfile import ZipFile

parser = argparse.ArgumentParser(description="Optimizes Kidney Exchange given by input file")
parser.add_argument('--inputFiles', nargs='+', default = ["data.dat"], help="list of .dat files to be used as input. List of number of \
                    incompatible pairs, number of compatible pairs, and list of all possible pairs with donor, recipient, egs")
parser.add_argument('--inputZipFile', help='filename of zip file to look for inputFiles in, if not given data is assumed to be uncompressed')
parser.add_argument('--quality', action='store_true', help="Optimize for quality")
parser.add_argument('-o', '--output')
parser.add_argument('--agents', help='output the quality of each agent to this file (.csv)')
parser.add_argument('--incompatibleOnly', action = 'store_true', help = "Run the oracle on only the incompatible pairs")
args=parser.parse_args()

data = []

def COUNT(v):
    if v[1] == 0:
        return 1
    if v[2] == 0:
        return 2
    return 3

results = ''
agentInfo = ''
if args.inputZipFile:
    inputZipFile = ZipFile(args.inputZipFile)
for fn in args.inputFiles:
    if args.inputZipFile:
        d = pickle.loads(inputZipFile.read(fn))
    else: 
        with open(fn, 'rb') as f:
            d = pickle.load(f)
    num_incompat = d[0]
    num_compat = d[1]
    num_pairs = num_incompat + num_compat
    matches = d[4]
    directed_matches = d[7]
    demo = d[5]
    departure_times = d[8]
    T = num_compat
    C = T
    model = pulp.LpProblem('match everyone', pulp.LpMaximize)
    if args.incompatibleOnly:
        matchVars = [v for v in matches if v[0] > T]
    else:
        matchVars = [v for v in matches]
    x = pulp.LpVariable.dicts('match',matchVars,lowBound=0, upBound=1, cat=pulp.LpInteger)
    matchVars = set(matchVars)
    if (args.quality):
        model += lpSum(x[v]*matches[v] for v in matchVars)

    else:
        model += lpSum(x[v]*COUNT(v) for v in matchVars)
    
    for t in range(1,num_pairs+1):
        model += lpSum(x[t,i,j] for i in range(num_incompat+1) for j in range(num_incompat+1) if (t,i,j) in matchVars) <= 1
    
    for i in range(1,num_incompat+1):
        model += lpSum(x[t,i,j] for t in range(1,num_pairs+1) for j in range(num_incompat+1) if (t,i,j) in matchVars) + \
                lpSum(x[t,j,i] for t in range(1,num_pairs+1) for j in range(1, num_incompat+1) if (t,j,i) in matchVars) + \
                lpSum(x[i+T,j,jp] for j in range(1,num_incompat+1) for jp in range(num_incompat+1) if (i+T,j,jp) in matchVars) <= 1, \
                'symetry ' + str(i)

        
    model.solve()
    
    if (args.incompatibleOnly):
        quality = sum(x[v].value()*matches[v] for v in matchVars) + sum(matches[(v,0,0)] for v in range(1,num_compat+1))
        count = sum(COUNT(v)*x[v].value() for v in matchVars) + num_compat
    else:
        quality = sum(x[v].value()*matches[v] for v in matchVars)
        count = sum(COUNT(v)*x[v].value() for v in matchVars)
    
    results += str(count) + "\t" + str(quality) + "\n"


    used_incompat = set()

    for v in matchVars:
        if round(x[v].value()) != 0:
            #if there is a compatible pair in the match
            if (v[0] <= T):
                #if compatible matched with itself
                if (v[1] == 0):
                    agentInfo += "C" + str(v[0]) + "\t" + str(0) + "\t" + str(directed_matches[v[0],0]) + "\t" \
                    + "C" + "\t" + str(directed_matches[v[0],0]) + "\t" + "C" + '\t' + str(demo[v[0]-1][20]) + '\t' + \
                    str(demo[v[0]-1][20]) + "\n"
                #compatible and incompatible
                elif (v[2] == 0):
                    agentInfo += "C" + str(v[0]) + "\t" + str(0) + "\t" + str(directed_matches[v[1]+T,v[0]]) + "\t" \
                    + "I" + "\t" + str(directed_matches[v[0],v[1]+T]) + "\t" + "I" + '\t' + str(demo[v[0]-1][20]) + '\t' +\
                    str(demo[v[0]-1][20]) + "\n"
                    agentInfo += "I" + str(v[1]) + "\t" + str(0) + "\t" + str(directed_matches[v[0],v[1]+T]) + "\t" \
                    + "C" + "\t" + str(directed_matches[v[1]+T,v[0]]) + "\t" + "C" + '\t' + str(demo[v[1]+C-1][20]) + '\t' + \
                    str(departure_times[v[1]-1]) + "\n"
                    used_incompat.add(v[1])
                #compatible and 2 incompatible
                else:
                    agentInfo += "C" + str(v[0]) + "\t" + str(0) + "\t" + str(directed_matches[v[2]+T,v[0]]) + "\t" \
                    + "I" + "\t" + str(directed_matches[v[0],v[1]+T]) + "\t" + "I" + '\t' + str(demo[v[0]-1][20]) + '\t' + \
                    str(demo[v[0]-1][20]) + "\n"
                    agentInfo += "I" + str(v[1]) + "\t" + str(0) + "\t" + str(directed_matches[v[0],v[1]+T]) + "\t" \
                    + "C" + "\t" + str(directed_matches[v[1]+T,v[2]+T]) + "\t" + "I" +  '\t' + str(demo[v[1]+C-1][20]) + '\t' +\
                    str(departure_times[v[1]-1]) + "\n"
                    agentInfo += "I" + str(v[2]) + "\t" + str(0) + "\t" + str(directed_matches[v[1]+T,v[2]+T]) + "\t" \
                    + "I" + "\t" + str(directed_matches[v[2]+T,v[0]]) + "\t" + "C" + '\t' + str(demo[v[2]+C-1][20]) + '\t' +\
                    str(departure_times[v[2]-1]) + "\n"
                    used_incompat.add(v[1])
                    used_incompat.add(v[2])
                  
            #only incompatible pairs
            else:
                #2 cycle
                if (v[2]==0):
                    agentInfo += "I" + str(v[0]-T) + "\t" + str(0) + "\t" + str(directed_matches[v[1]+T,v[0]]) + "\t" \
                    + "I" + "\t" + str(directed_matches[v[0],v[1]+T]) + "\t" + "I" + '\t' + str(demo[v[0]-1][20]) + '\t' +\
                    str(departure_times[v[0]-C-1]) + "\n"
                    agentInfo += "I" + str(v[1]) + "\t" + str(0) + "\t" + str(directed_matches[v[0],v[1]+T]) + "\t" \
                    + "I" + "\t" + str(directed_matches[v[1]+T,v[0]]) + "\t" + "I" + '\t' + str(demo[v[1]-1][20]) + '\t' + \
                    str(departure_times[v[1]-1]) + "\n"  
                    used_incompat.add(v[0]-T)
                    used_incompat.add(v[1])
                #3 cycle
                else:
                    agentInfo += "I" + str(v[0]-T) + "\t" + str(0) + "\t" + str(directed_matches[v[2]+T,v[0]]) + "\t" \
                    + "I" + "\t" + str(directed_matches[v[0],v[1]+T]) + "\t" + "I" + '\t' + str(demo[v[0]-1][20]) + '\t' +\
                    str(departure_times[v[0]-C-1]) + "\n"
                    agentInfo += "I" + str(v[1]) + "\t" + str(0) + "\t" + str(directed_matches[v[0],v[1]+T]) + "\t" \
                    + "I" + "\t" + str(directed_matches[v[1]+T,v[2]+T]) + "\t" + "I" + '\t' + str(demo[v[1]+C-1][20]) + '\t'+\
                    str(departure_times[v[1]-1]) + "\n"
                    agentInfo += "I" + str(v[2]) + "\t" + str(0) + "\t" + str(directed_matches[v[1]+T,v[2]+T]) + "\t" \
                    + "I" + "\t" + str(directed_matches[v[2]+T,v[0]]) + "\t" + "I" + '\t' + str(demo[v[2]+C-1][20]) + '\t' + \
                    str(departure_times[v[2]-1]) + "\n"
                    used_incompat.add(v[0]-T)
                    used_incompat.add(v[1])
                    used_incompat.add(v[2])

    for i in range(1,num_incompat+1):
        if i not in used_incompat:
             agentInfo += "I" + str(i) + "\t" + str(d[2]+2) + "\t" + str(0) + "\t" \
        + "N" + "\t" + str(0) + "\t" + "N" + '\t' + str(demo[i+C-1][20]) + '\t' + str(departure_times[i-1]) + "\n"

if args.inputZipFile:
    inputZipFile.close()

if args.output:
    with open(args.output, 'w') as f:
        f.write(results)
else:
    print results

if args.agents:
    with open(args.agents, 'w') as f:
        f.write(agentInfo)
