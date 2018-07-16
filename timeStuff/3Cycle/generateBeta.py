from gurobipy import *
import json
import pickle
import argparse
import os

parser = argparse.ArgumentParser(description="Computes the dual of our problem")
parser.add_argument('--inputFiles', nargs='+', default = ["data.dat"], help='input file to use')
parser.add_argument('-o', '--output', help='output file (json) to use. Includes demographic information and beta values')
parser.add_argument('--quality', action = "store_true", help="Optimize for quality")

args = parser.parse_args()

def COUNT(v):
    if v[2] != 0:
        return 3
    if v[1] != 0:
        return 2
    return 1

results = []
    
for fn in args.inputFiles:
    with open(fn, 'rb') as f:
        d = pickle.load(f)
    #T is number of compatible pairs
    T = d[1]
    #K is number of incompatible pairs
    K = d[0]
    matches = d[4]
    demo = d[5]
    
    model = Model("Dual Optimizer")
    
    alpha = {}
    beta = {}
    
    for t in range(1,T+1):
        if any(k[0] == t for k in matches):
            alpha[t] = model.addVar(vtype = GRB.CONTINUOUS, lb=0, name='alpha_'+str(t))
    
    for i in range(1,K+1):
        if any(k[0] == i+T for k in matches):
            beta[i] = model.addVar(vtype=GRB.CONTINUOUS, lb=0, name='beta_'+str(i))
        elif any(k[1] == i for k in matches):
            beta[i] = model.addVar(vtype=GRB.CONTINUOUS, lb=0, name='beta_'+str(i))
        elif any(k[2] == i for k in matches):
            beta[i] = model.addVar(vtype=GRB.CONTINUOUS, lb=0, name='beta_'+str(i))

    beta[0] = 0
    
    if args.quality:
        model.addConstrs((matches[t,i,j] - (alpha[t] if t in alpha else 0) - beta[i] - beta[j] - (beta[t-T] if t-T in beta else 0) <= 0 for t in range(1,T+K+1) if t in alpha or t-T in beta for i in beta  \
                 for j in beta if (t,i,j) in matches), "something...")
    else:
        model.addConstrs((COUNT((t,i,j)) - (alpha[t] if t in alpha else 0) - beta[i] - beta[j] - (beta[t-T] if t-T in beta else 0) <= 0 for t in range(1, T+K+1) if t in alpha or t-T in beta for i in beta \
                for j in beta if (t,i,j) in matches), "something...")
    
    obj = quicksum(alpha[t] for t in alpha) + quicksum(beta[i] for i in beta)
    model.setObjective(obj, GRB.MINIMIZE)
    model.optimize()
    
    for i in range(1,K+1):
        results.append((demo[i+T-1], (beta[i].X if i in beta else 0)))

if args.output:
    with open(args.output, 'w') as f:
        f.write(json.dumps(results))
else:
    print json.dumps(results)
