#generate data that can be used to create 2 cycle loops (or infinite size loops or incompatible pairs)
import random
import json
import argparse

parser = argparse.ArgumentParser(description="Generates data for compatible and incompatible kidney exchange pairs")
parser.add_argument('-o', '--outputFile', nargs = '?', default = "data.json", help = "JSON File for data to be saved to")
parser.add_argument('-T', '--numAgents', nargs = '?', default = 100, help = "Number of agents", dest="T", type=int)
parser.add_argument('-m', '--mean', nargs = '?', default = 20, help = "Mean quality", dest="mean", type=int)
args = parser.parse_args()

T = args.T
mean = args.mean
#generate list of compatible pairs
compat = []
#d=0 is more compatible and d=1 is less compatible, r=0 is easy to match, r=1 hard to match
for i in range(T):
    if random.random() < .6:
        d = 0
    else:
        d = 1
    if random.random() < .5:
        r = 1
    else:
        r = 0
    q = random.expovariate(1./mean)
    compat.append((d,r,q))

#generate list of incompatible pairs  
incompat = []        
for i in range(T):
    if random.random() < .5:
        d = 0
    else:
        d = 1
    if random.random() < .3:
        r = 1
    else:
        r = 0
    q = random.expovariate(1./mean)
    incompat.append((d,r,q))

#create matrix with 1 representing an edge between the incompatible pairs
Igraph = []
for i in range(T):
    Igraph.append([])
    for j in range(T):
        if i==j:
            Igraph[i].append(0)
            continue
        p = random.random()
        if not incompat[i][0] and not incompat[j][1]:
            Igraph[i].append(int(p < .04))
        if not incompat[i][0] and incompat[j][1]:
            Igraph[i].append(int(p < .03))
        if incompat[i][0] and not incompat[j][1]:
            Igraph[i].append(int(p < .02))
        if incompat[i][0] and incompat[j][1]:
            Igraph[i].append(int(p < .01))

#create matrix with 1 representing an edge between compatible and incompatible pair
CIgraph = []
for i in range(T):
    CIgraph.append([])
    for j in range(T):
        p = random.random()
        p2 = random.random()
        v1 = 0
        v2 = 0
        if not compat[i][0] and not incompat[j][1]:
            v1 = int(p < .04)
        if not compat[i][0] and incompat[j][1]:
            v1 = int(p < .03)
        if compat[i][0] and not incompat[j][1]:
            v1 = int(p < .02)
        if compat[i][0] and incompat[j][1]:
            v1 = int(p < .01)
        #add v2 because the graph is undirected
        if not incompat[j][0] and not compat[i][1]:
            v2 = int(p2 < .04)
        if not incompat[j][0] and compat[i][1]:
            v2 = int(p2 < .03)
        if incompat[j][0] and not compat[i][1]:
            v2 = int(p2 < .02)
        if incompat[j][0] and compat[i][1]:
            v2 = int(p2 < .01)
        CIgraph[i].append(v1&v2)
            

with open(args.outputFile, 'w') as f:
    f.write(json.dumps((compat, incompat, Igraph, CIgraph)))
