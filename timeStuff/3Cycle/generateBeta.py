from gurobipy import *
import json
import pickle
import argparse
import os

parser = argparse.ArgumentParser(description="Computes the dual of our problem")
parser.add_argument('--inputFiles', nargs='+', default = ["data.dat"], help='input file to use')
parser.add_argument('-o', '--output', help='output file (json) to use. Includes demographic information and beta values')
parser.add_argument('--quality', action = "store_true", help="Optimize for quality")
parser.add_argument('--graph_state', action='store_true', help='Flag should be present if online LP estimation is included in training data')

args = parser.parse_args()

def COUNT(v):
    if v[2] != 0:
        return 3
    if v[1] != 0:
        return 2
    return 1

def calcBetaLP(C, matches, available_incompat):
    estimator = Model("estimate beta vals")
    beta = {}
    for i in available_incompat:
        if any(k[0] == i+C for k in matches if k[1] in available_incompat and k[2]  in available_incompat):
            beta[i] = estimator.addVar(vtype=GRB.CONTINUOUS, lb=0, name='beta_'+str(i))
        elif any(k[1] == i for k in matches if k[0] > C and k[0]-C  in available_incompat and k[2] in available_incompat):
            beta[i] = estimator.addVar(vtype=GRB.CONTINUOUS, lb=0, name='beta_'+str(i))
        elif any(k[2] == i for k in matches if k[0] > C and k[0]-C  in available_incompat and k[1]  in available_incompat):
            beta[i] = estimator.addVar(vtype=GRB.CONTINUOUS, lb=0, name='beta_'+str(i))
    if args.quality:
        estimator.addConstrs((matches[t+C,i,j] - beta[i] -beta[j] - beta[t]  <= 0 for t in beta  for i in beta \
                for j in beta if (t+C,i,j) in matches), 'something')
    else:
        estimator.addConstrs((COUNT((t+C,i,j)) - beta[i] - beta[j] - beta[t] <= 0 for t in beta  for i in beta \
                for j in beta if (t+C,i,j) in matches), 'something')
    obj = quicksum(beta[i] for i in beta)
    estimator.setObjective(obj, GRB.MINIMIZE)
    estimator.optimize()
    newBeta = {i:beta[i].X for i in beta}
    for i in available_incompat:
        if i not in newBeta:
            newBeta[i] = 0
    return newBeta


results = []
    
for fn in args.inputFiles:
    with open(fn, 'rb') as f:
        d = pickle.load(f)
    #T is number of compatible pairs
    C = d[1]
    #K is number of incompatible pairs
    K = d[0]
    T = d[2]
    matches = d[4]
    demo = d[5]
    departure_times = d[8]
    
    model = Model("Dual Optimizer")
    
    alpha = {}
    beta = {}
    
    for t in range(1,C+1):
        if any(k[0] == t for k in matches):
            alpha[t] = model.addVar(vtype = GRB.CONTINUOUS, lb=0, name='alpha_'+str(t))
    
    for i in range(1,K+1):
        if any(k[0] == i+C for k in matches):
            beta[i] = model.addVar(vtype=GRB.CONTINUOUS, lb=0, name='beta_'+str(i))
        elif any(k[1] == i for k in matches):
            beta[i] = model.addVar(vtype=GRB.CONTINUOUS, lb=0, name='beta_'+str(i))
        elif any(k[2] == i for k in matches):
            beta[i] = model.addVar(vtype=GRB.CONTINUOUS, lb=0, name='beta_'+str(i))

    beta[0] = 0
    
    if args.quality:
        model.addConstrs((matches[t,i,j] - (alpha[t] if t in alpha else 0) - beta[i] - beta[j] - (beta[t-C] if t-C in beta else 0) <= 0 for t in range(1,C+K+1) if t in alpha or t-C in beta for i in beta  \
                 for j in beta if (t,i,j) in matches), "something...")
    else:
        model.addConstrs((COUNT((t,i,j)) - (alpha[t] if t in alpha else 0) - beta[i] - beta[j] - (beta[t-C] if t-C in beta else 0) <= 0 for t in range(1, C+K+1) if t in alpha or t-C in beta for i in beta \
                for j in beta if (t,i,j) in matches), "something...")
    
    obj = quicksum(alpha[t] for t in alpha) + quicksum(beta[i] for i in beta)
    model.setObjective(obj, GRB.MINIMIZE)
    model.optimize()
    
    if args.graph_state:
        available_incompat = set()
        arriving_incompat = {}
        arriving_compat = {}


        for i in range(C):
            if demo[i][20] not in arriving_compat:
                arriving_compat[demo[i][20]] = []
            arriving_compat[demo[i][20]].append(i+1)
        for i in range(C,C+K):
            if demo[i][20] not in arriving_incompat:
                arriving_incompat[demo[i][20]] = set()
            arriving_incompat[demo[i][20]].add(i+1-C)

        for t in range(T):
            departing_incompat = set()
            for i in available_incompat:
                if departure_times[i-1] < t:
                    departing_incompat.add(i)
            available_incompat = available_incompat.difference(departing_incompat)

            if t in arriving_incompat:
                available_incompat = available_incompat.union(arriving_incompat[t])
            lpBeta = calcBetaLP(C, matches, available_incompat)
            for i in available_incompat:
                results.append((demo[i+C-1], (beta[i].X if i in beta else 0), lpBeta[i]))
    else:
        for i in range(1,K+1):
            results.append((demo[i+C-1], (beta[i].X if i in beta else 0)))

if args.output:
    with open(args.output, 'w') as f:
        f.write(json.dumps(results))
else:
    print json.dumps(results)
