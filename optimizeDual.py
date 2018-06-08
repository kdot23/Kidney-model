import json
from gurobipy import *
import numpy as np

with open("data.json",'r') as f:
    data = json.load(f)
compat = data[0]
incompat = data[1]
Igraph = data[2]
CIgraph = data[3]
T = len(compat)

model = Model('Dual Kidney Optimizer')
alpha = {}
beta = {}
gamma = {}
for i in range(T):
    alpha[i] = model.addVar(vtype=GRB.CONTINUOUS, name="alpha"+str(i))
for i in range(T):
    if sum(Igraph[i]) or sum(1 for j in range(T) if CIgraph[j][i] and incompat[i][2]>compat[j][2]) > 0:
        beta[i] = model.addVar(vtype=GRB.CONTINUOUS, lb=0, name="beta"+str(i))
for i in range(T):
    for j in range(T):
        if Igraph[i][j] and Igraph[j][i]:
            gamma[(i,j)] =  model.addVar(vtype=GRB.CONTINUOUS, name="gamma"+str(i) + ", " + str(j))

model.addConstrs((alpha[i]+beta[j]-compat[i][2]-incompat[j][2] >= 0 for i in alpha for j in beta), name = "something???")
model.addConstrs((alpha[i]-compat[i][2] >= 0 for i in alpha), name = "something2")
model.addConstrs((beta[v[0]] + gamma[v] - gamma[(v[1],v[0])] - incompat[v[0]][2] >= 0 for v in gamma ), name = "something3")

obj = quicksum(alpha[i] for i in alpha) + quicksum(beta[i] for i in beta)
model.setObjective(obj, GRB.MINIMIZE)
model.optimize()

for v in model.getVars():
    if v.X != 0:
        print str(v.Varname) + ' ' + str(v.X)
print obj.getValue() 