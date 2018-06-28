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

results = []

data = []
for fn in args.inputFiles:
    with open (fn, 'rb') as f:
        data.append(pickle.load(f))

for d in data:
    #T is number of compatible pairs
    T = d[0]
    #K is number of incompatible pairs
    K = d[1]
    matches = d[2]
    demo = d[4]
    
    model = Model("Dual Optimizer")
    
    alpha = {}
    beta = {}
    
    for t in range(1,T+K+1):
        if any(k[0] == t for k in matches):
            alpha[t] = model.addVar(vtype = GRB.CONTINUOUS, lb=0, name='alpha_'+str(t))
    
    for i in range(1,K+1):
        if any(k[0] == i+T for k in matches):
            beta[i] = model.addVar(vtype=GRB.CONTINUOUS, lb=0, name='beta_'+str(i))
        if any(k[1] == i for k in matches):
            beta[i] = model.addVar(vtype=GRB.CONTINUOUS, lb=0, name='beta_'+str(i))
    beta[0] = 0
    
    if args.quality:
        model.addConstrs((matches[t,i] - alpha[t] - beta[i] - (beta[t-T] if t-T in beta else 0) <= 0 for t in alpha for i in beta  \
                if (t,i) in matches), "something...")
    else:
        model.addConstrs((1 - alpha[t] - beta[i] - (beta[t-T] if t-T in beta else 0) <= 0 for t in alpha for i in beta \
                if (t,i) in matches), "something...")
    
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
