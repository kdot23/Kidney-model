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
iMatches = {}
ciMatches = {}
#only create variables for edges which exist
for i in range(T):
    for j in range(T):
        if Igraph[i][j] and Igraph[j][i]:
            iMatches[(i,j)] = model.addVar(vtype = GRB.CONTINUOUS, lb = 0, ub = 1, name = "I" + str((i,j)))
        if CIgraph[i][j] and incompat[j][2] > compat[i][2]:
            ciMatches[(i,j)] = model.addVar(vtype = GRB.CONTINUOUS, lb = 0, ub = 1, name = "CI" + str((i,j)))
    ciMatches[(i,-1)] = model.addVar(vtype = GRB.CONTINUOUS, lb = 0, ub = 1, name = "CI" + str((i,-1)))

#first entry is which compatible pair, second entry is incompatible pair or -1 for self matching
#ciMatches are compatible to incompatible or compatible to themselves (-1)
#ciMatches = model.addVars(range(T), range(-1, T), vtype = GRB.BINARY, name="CIMatches")
#iMatches are incompatible to incompatible
#iMatches = model.addVars(range(T), range(T), vtype = GRB.BINARY, name = "IMatches")

#compatible pair can only match with up to one incompatible pair or themselves
model.addConstrs((quicksum(ciMatches[i,j] for j in range(-1,T) if (i,j) in ciMatches)  == 1 for i in range(T)), "Compat to Incompat Matches")
#incompatible pair can only match with up to one incompatible pair or one compatible pair
model.addConstrs((quicksum(iMatches[i,j] for j in range(T) if (i,j) in iMatches) + quicksum(ciMatches[j,i] for j in range(T) if (j,i) in ciMatches) <= 1 for i in range(T)), "Incompat Matches")

model.addConstrs((iMatches[i,j] == iMatches[j,i] for i in range(T) for j in range(T) if (i,j) in iMatches), "Symmetry")

#compatible will only exchange if they have a higher quality score due to the exchange
#model.addConstrs((ciMatches[i,j]*incompat[j][2] >= ciMatches[i,j]*compat[i][2] for i in range(T) for j in range(T)), " Compat Selfishness")

#an edge between C and I can only exist if the pairs are compatible
#model.addConstrs((CIgraph[i][j] >= ciMatches[i,j] for i in range(T) for j in range(T)), "Compatible to Incompatible Compatibility")

#an edge between two Is can only exist if the pairs are compatible 
#model.addConstrs((Igraph[i][j] >= iMatches[i,j] for i in range(T) for j in range(T)), "Incompatible Compatibility")

#substitute this line for line 24 to allow for arbitrary length cycles in incompatible pool
#model.addConstrs((quicksum(iMatches[i,j] for j in range(T) if (i,j) in iMatches) == quicksum(iMatches[j,i] for j in range(T)) for i in range(T)), "Give one to get one")

#obj = quicksum(ciMatches[i,j] for i in range(T) for j in range(T) if (i,j) in ciMatches) + quicksum(ciMatches[i,-1] for i in range(T) ) + quicksum(iMatches[i,j] for i in range(T) for j in range(T) if (i,j) in iMatches)
obj = quicksum((compat[i][2]+incompat[j][2])*ciMatches[i,j] for i in range(T) for j in range(T) if (i,j) in ciMatches) + quicksum(compat[i][2]*ciMatches[i,-1] for i in range(T) ) + quicksum(incompat[i][2]*iMatches[i,j] for i in range(T) for j in range(T) if (i,j) in iMatches)
model.setObjective(obj, GRB.MAXIMIZE) 

model.optimize()
count = 0
for v in model.getVars():
    if v.X != 0:
        print str(v.Varname) + ' ' + str(v.X)
        count+=1
print count
util = 0.
count2 = 0
count3 = 0
count4 = 0
for var in iMatches:
    if (iMatches[var].X != 0):
        util += incompat[var[0]][2]
        count4 += 1
for var in ciMatches:
    if (ciMatches[var].X != 0):
        util += compat[var[0]][2]+incompat[var[1]][2]
        count3+=1
        if var[1] != -1:
            count2+=1
            count4+=1
        
print util
print float(count2)/count3
print float(count4)
print(obj.getValue())