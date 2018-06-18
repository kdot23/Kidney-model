#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 15 11:05:58 2018

@author: kelseylieberman
"""
import json
import argparse
from gurobipy import *

parser = argparse.ArgumentParser(description="Optimizes Kidney Exchange given by input file")
parser.add_argument('--inputFile', nargs='?', help="JSON File to be used as input. List of number of \
                    incompatible pairs, number of compatible pairs, and list of all possible pairs with donor, recipient, egs")
parser.add_argument('--quality', action='store_true', help="Optimize for quality")
args=parser.parse_args()
with open (args.inputFile, 'r') as f:
    data = json.load(f)


num_incompat = data[0]
num_compat = data[1]
num_pairs = num_incompat + num_compat
matches = data[2]
T = num_compat
model = Model('Kideny Optimizer')
matchVars = {}
for i in range(num_pairs):
    for j in range(num_incompat+1):
        if matches[i][j] != 0:
            matchVars[(i,j)] = model.addVar(vtype = GRB.CONTINUOUS, lb = 0, ub=1,  name = "match_" + str((i,j)))

model.addConstrs((quicksum(matchVars[t,i] for i in range(num_incompat+1) if (t,i) in matchVars) <= 1 \
                  for t in range(num_pairs)), "Only match with one pair")

model.addConstrs((quicksum(matchVars[t,i] for t in range(num_pairs) if (t,i) in matchVars) + \
                  quicksum(matchVars[i+T-1,j] for j in range(1,num_incompat+1) if (i+T-1,j) in matchVars) <= 1 for i in range(1,num_incompat+1)), \
                    "undirected graph")
#add symmetry constraint??
if (args.quality):
    obj = quicksum(matchVars[i,j]*matches[i][j] for i in range(num_pairs) for j in range(num_incompat+1) if (i,j) in matchVars)
else:
    obj = quicksum(matchVars[i,j] for i in range(num_pairs) for j in range(num_incompat+1) if (i,j) in matchVars)
    
model.setObjective(obj, GRB.MAXIMIZE) 
model.optimize()
