#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 15 11:05:58 2018

@author: kelseylieberman
"""

parser = argparse.ArgumentParser(description="Optimizes Kidney Exchange given by input file")
parser.add_argument('--inputFile', nargs='?', help="JSON File to be used as input. List of number of \
                    incompatible pairs, number of compatible pairs, and list of all possible pairs with donor, recipient, egs")
parser.add_argument('--quality', action='store_true', help="Optimize for quality")
args=parser.parse_args()
num_incompat = args.inputFile[0]
num_compat = args.inputFile[1]
num_pairs = num_incompat + num_compat
matches = args.inputFile[2]

model = Model('Kideny Optimizer')
matchVars = {}
for i in range(num_pairs):
    for j in range(num_pairs):
        if matches[i][j][2] != 0:
            matchVars[(i,j)] = model.addVar(vtype = GRB.CONTINUOUS, lb = 0, ub = 1, name = "match_" + str((i,j)))

model.addConstrs((quicksum(matchVars[i,j] for j in range(num_pairs) if (i,j) in matchVars) <= 1 \
                  for i in range(num_pairs)), "Compat  Matches")

if (args.quality):
    obj = quicksum(matchVars[i,j]*matches[i][j][2] for i in num_pairs for j in num_pairs)
else:
    obj = quicksum(matchVars[i,j] for i in num_pairs)
    
model.setObjective(obj, GRB.MAXIMIZE) 
model.optimize()