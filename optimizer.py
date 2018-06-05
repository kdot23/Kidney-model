import json
from gurobipy import *

with open("data.json",'r') as f:
    data = json.load(f)
compat = data[0]
imcompat = data[1]
Igraph = data[2]
Cgraph = data[3]
T = len(compat)

model = Model('Kideny Optimizer')
#first entry is which compatible pair, second entry is incompatible pair or -1 for self matching
#cMatches are compatible to incompatible
cMatches = model.addVars(range(T), range(-1, T), vtype = GRB.BINARY, name="CMatches")
#iMatches are incompatible to incompatible
iMatches = model.addVars(range(T), range(T), vtype = GRB.BINARY, name = "IMatches")

#compatible pair can only match with one incompatible pair
model.addConstrs((quicksum(cMatches[i,j] for j in range(-1,T))  <= 1 for i in range(T)), "Compat Matches")
#incompatible pair can only match with either one incompatible pair or one compatible pair
model.addConstrs((quicksum(iMatches[i,j] for j in range(T)) + quicksum(cMatches[j,i] for j in range(T)) <= 1 for i in range(T)), "Imcompat Matches")

model.addConstrs((iMatches[i,j] == iMatches[j,i] for i in range(T) for j in range(T)), "Symetry")

model.addConstrs((cMatches[i,j]*imcompat[j][2] >= cMatches[i,j]*compat[i][2] for i in range(T) for j in range(T)), " Compat Selfishness")

model.addConstrs((Cgraph[i][j] >= cMatches[i,j] for i in range(T) for j in range(T)), "Compatible Compatibility")

model.addConstrs((Igraph[i][j] >= iMatches[i,j] for i in range(T) for j in range(T)), "Incompatible Compatibility")

obj = quicksum((compat[i][2]+imcompat[j][2])*cMatches[i,j] for i in range(T) for j in range(T)) + quicksum(compat[i][2]*cMatches[i,-1] for i in range(T)) + quicksum(imcompat[i][2]*iMatches[i,j] for i in range(T) for j in range(T))
model.setObjective(obj, GRB.MAXIMIZE) 

model.optimize()
for v in model.getVars():
    if v.X != 0:
        print str(v.Varname) + ' ' + str(v.X)