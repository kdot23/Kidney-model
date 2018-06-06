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
CIgraph = data[3]
ICgraph = data[4]
T = len(compat)

def twoCycle(Igraph, CIgraph, ICgraph):
    twoCycles = {}
    for i in range(T):
        for j in range(T):
            #if an edge exists between incompatible pair i and j 
            if Igraph[i][j] and Igraph[j][i]:
                #value at (i,j) is the quality of both exchanges
                twoCycles[(i,j)] = incompat[i][2] + incompat[j][2]
            if CIgraph[i][j] and ICgraph[j][i] and incompat[j][2] > compat[i][2]:
                twoCycles[(i+T,j)] = compat[i][2] + incompat[j][2]
        twoCycles[(i+T,i+T)] = compat[i][2]
    return twoCycles
def threeCycle(Igraph, CIgraph, ICgraph):
    threeCycles = {}
    for i in range(T):
        for j in range(T):
            if Igraph[i][j]:
                for w in range(T):
                    if w >= i or w >= j:
                        break
                    if Igraph[j][w] and Igraph[w][i]:
                        threeCycles[(i,j,w)] = incompat[i][2] + incompat[j][2] + incompat[w][2]
            if CIgraph[i][j]:
                for w in range(T):
                    if w >= j:
                        break
                    if Igraph[j][w] and ICgraph[w][i]:
                        threeCycles[(i+T,j,w)] = compat[i][2] + incompat[j][2] + incompat[w][2]
    return threeCycles
                        

twoCycles = twoCycle(Igraph, CIgraph, ICgraph)
threeCycles = threeCycle(Igraph, CIgraph, ICgraph)
    



model = Model('Kideny Optimizer')
model2 = Model('Kidney Optimizer2')
c= {}
for cycle in twoCycles:
    c[cycle] = model2.addVar(vtype=GRB.BINARY, name="c_%s" % str(cycle))
for cycle in threeCycles:
    c[cycle] = model2.addVar(vtype=GRB.BINARY, name="c_%s" % str(cycle))
model2.update()

for v in range(2*T):
  constraint = []
  for cycle in c:
      if (v in cycle):
          constraint.append(c[cycle])
  if constraint:
      model2.addConstr( quicksum( con for con in constraint ) <= 1 , name="v%d" % v)
model2.setObjective( quicksum( c[cycle] * twoCycles[cycle] for cycle in twoCycles ) + quicksum(c[cycle] * threeCycles[cycle] for cycle in threeCycles ), GRB.MAXIMIZE )


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

model.addConstrs((CIgraph[i][j] >= cMatches[i,j] for i in range(T) for j in range(T)), "Compatible Compatibility")

model.addConstrs((Igraph[i][j] >= iMatches[i,j] for i in range(T) for j in range(T)), "Incompatible Compatibility")

model.addConstrs((ICgraph[i][j] >= icMatches[i,j] for i in range(T) for j in range(T)), "Incompatible to compatible compatibility")

obj = quicksum(compat[i][2]*cMatches[i,j] for i in range(T) for j in range(T)) + quicksum(incompat[j][2]*icMatches[j,i] for i in range(T) for j in range(T)) + quicksum(compat[i][2]*cMatches[i,-1] for i in range(T)) + quicksum(incompat[i][2]*iMatches[i,j] for i in range(T) for j in range(T))
model.setObjective(obj, GRB.MAXIMIZE) 

model2.optimize()
for v in model2.getVars():
    if v.X != 0:
        print str(v.Varname) + ' ' + str(v.X)