#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  5 14:30:40 2018

@author: kelseylieberman
"""

import json
from gurobipy import *

with open("data.json",'r') as f:
    data = json.load(f)
compat = data[0]
incompat = data[1]
Igraph = data[2]
Cgraph = data[3]
ICgraph = data[4]
T = len(compat)

model = Model('Kideny Optimizer')
#first entry is which compatible pair, second entry is incompatible pair or -1 for self matching
#cMatches are compatible to incompatible
cMatches = model.addVars(range(T), range(-1, T), vtype = GRB.BINARY, name="CMatches")
#iMatches are incompatible to incompatible
iMatches = model.addVars(range(T), range(T), vtype = GRB.BINARY, name = "IMatches")
icMatches = model.addVars(range(T), range(T), vtype = GRB.BINARY, name = "ICMatches")

#compatible pair can only match with one incompatible pair
model.addConstrs((quicksum(cMatches[i,j] for j in range(-1,T))  <= 1 for i in range(T)), "Compat Matches")
#incompatible pair can only match with either one incompatible pair or one compatible pair
model.addConstrs((quicksum(iMatches[i,j] for j in range(T)) + quicksum(icMatches[i,j] for j in range(T)) <= 1 for i in range(T)), "Incompat Matches")

#model.addConstrs((iMatches[i,j] == iMatches[j,i] for i in range(T) for j in range(T)), "Symmetry")
model.addConstrs((quicksum(iMatches[i,j] for j in range(T)) + quicksum(icMatches[i,j] for j in range(T)) == quicksum(iMatches[j,i] for j in range(T)) + quicksum(cMatches[j,i] for j in range(T)) for i in range(T)), "Give one to get one")

model.addConstrs((quicksum(cMatches[i,j] for j in range(T)) == quicksum(icMatches[j,i] for j in range(T)) for i in range(T)), "Give one to get one compat")


model.addConstrs((cMatches[i,j]*incompat[j][2] >= cMatches[i,j]*compat[i][2] for i in range(T) for j in range(T)), " Compat Selfishness")

model.addConstrs((Cgraph[i][j] >= cMatches[i,j] for i in range(T) for j in range(T)), "Compatible Compatibility")

model.addConstrs((Igraph[i][j] >= iMatches[i,j] for i in range(T) for j in range(T)), "Incompatible Compatibility")

model.addConstrs((ICgraph[i][j] >= icMatches[i,j] for i in range(T) for j in range(T)), "Incompatible to compatible compatibility")

obj = quicksum(compat[i][2]*cMatches[i,j] for i in range(T) for j in range(T)) + quicksum(incompat[j][2]*icMatches[j,i] for i in range(T) for j in range(T)) + quicksum(compat[i][2]*cMatches[i,-1] for i in range(T)) + quicksum(incompat[i][2]*iMatches[i,j] for i in range(T) for j in range(T))
model.setObjective(obj, GRB.MAXIMIZE) 

model.optimize()
for v in model.getVars():
    if v.X != 0:
        print str(v.Varname) + ' ' + str(v.X)