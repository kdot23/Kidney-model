import json
from gurobipy import *

with open("data.json",'r') as f:
    data = json.load(f)
compat = data[0]
incompat = data[1]
Igraph = data[2]
CIgraph = data[3]
T = len(compat)

model = Model('Kideny Optimizer')
#first entry is which compatible pair, second entry is incompatible pair or -1 for self matching
#ciMatches are compatible to incompatible or compatible to themselves (-1)
ciMatches = model.addVars(range(T), range(-1, T), vtype = GRB.BINARY, name="CIMatches")
#iMatches are incompatible to incompatible
iMatches = model.addVars(range(T), range(T), vtype = GRB.BINARY, name = "IMatches")

#compatible pair can only match with up to one incompatible pair or themselves
model.addConstrs((quicksum(ciMatches[i,j] for j in range(-1,T))  <= 1 for i in range(T)), "Compat to Incompat Matches")
#incompatible pair can only match with up to one incompatible pair or one compatible pair
model.addConstrs((quicksum(iMatches[i,j] for j in range(T)) + quicksum(ciMatches[j,i] for j in range(T)) <= 1 for i in range(T)), "Incompat Matches")

#model.addConstrs((iMatches[i,j] == iMatches[j,i] for i in range(T) for j in range(T)), "Symmetry")

#compatible will only exchange if they have a higher quality score due to the exchange
model.addConstrs((ciMatches[i,j]*incompat[j][2] >= ciMatches[i,j]*compat[i][2] for i in range(T) for j in range(T)), " Compat Selfishness")

#an edge between C and I can only exist if the pairs are compatible
model.addConstrs((CIgraph[i][j] >= ciMatches[i,j] for i in range(T) for j in range(T)), "Compatible to Incompatible Compatibility")

#an edge between two Is can only exist if the pairs are compatible 
model.addConstrs((Igraph[i][j] >= iMatches[i,j] for i in range(T) for j in range(T)), "Incompatible Compatibility")

#substitute this line for line 24 to allow for arbitrary length cycles in incompatible pool
model.addConstrs((quicksum(iMatches[i,j] for j in range(T)) == quicksum(iMatches[j,i] for j in range(T)) for i in range(T)), "Give one to get one")

obj = quicksum((compat[i][2]+incompat[j][2])*ciMatches[i,j] for i in range(T) for j in range(T)) + quicksum(compat[i][2]*ciMatches[i,-1] for i in range(T)) + quicksum(incompat[i][2]*iMatches[i,j] for i in range(T) for j in range(T))
model.setObjective(obj, GRB.MAXIMIZE) 

model.optimize()
for v in model.getVars():
    if v.X != 0:
        print str(v.Varname) + ' ' + str(v.X)