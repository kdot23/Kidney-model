from gurobipy import *
import json
import pickle
import argparse
import os

parser = argparse.ArgumentParser(description="Computes the dual of our problem")
parser.add_argument('--inputFiles', nargs='+', default = ["data.dat"], help='input file to use')
parser.add_argument('-o', '--output', help='output file (json) to use. Includes demographic information and beta values')
parser.add_argument('--quality', action = "store_true", help="Optimize for quality")
parser.add_argument('--graph_state', action='store_true', help='present if onlineLP estimation of beta should be included in training data')

args = parser.parse_args()

results = []

def COUNT(v):
    if v[1] == 0:
        return 1
    return 2
def calcBetasLP(T, K, matches):
    estimator = Model('estimate beta values')
    alpha = {}
    beta = {}
    for i in range(1,K+1):
        if any(k[0] == i+T for k in matches) or any(k[1] == i for k in matches if k[0] > T ):
            beta[i] = estimator.addVar(vtype=GRB.CONTINUOUS, lb=0, \
                    name='beta_'+str(i))
    if args.quality:
        estimator.addConstrs((matches[t,i] -  beta[i] -  (beta[t-T] if t-T in beta else 0) <= 0 \
                for t in range(T+1,T+K+1) if t-T in beta for i in beta if (t,i) in matches),  'something...')
    else:
        estimator.addConstrs((COUNT((t,i)) - beta[i] - (beta[t-T] if t-T in beta else 0) <= 0 \
                for t in range(T+1, T+K+1) if t-T in beta for i in beta if (t,i) in matches), 'something...')
    obj = quicksum(beta[i] for i in beta)
    estimator.setObjective(obj, GRB.MINIMIZE)
    estimator.optimize()
    newBeta =  {i:beta[i].X for i in beta}
    for i in range(1, K+1):
        if i not in newBeta:
            newBeta[i] = 0
    return newBeta

for fn in args.inputFiles:
    with open(fn, 'rb') as f:
        d = pickle.load(f)
    #T is number of compatible pairs
    K = d[0]
    #K is number of incompatible pairs
    T = d[1]
    matches = d[2]
    demo = d[4]
    if args.graph_state:
        lpbeta = calcBetasLP(T,K,matches)
    
    model = Model("Dual Optimizer")
    
    alpha = {}
    beta = {}
    
    for t in range(1,T+1):
        if any(k[0] == t for k in matches):
            alpha[t] = model.addVar(vtype = GRB.CONTINUOUS, lb=0, name='alpha_'+str(t))
    
    for i in range(1,K+1):
        if any(k[0] == i+T for k in matches):
            beta[i] = model.addVar(vtype=GRB.CONTINUOUS, lb=0, name='beta_'+str(i))
            continue
        if any(k[1] == i for k in matches):
            beta[i] = model.addVar(vtype=GRB.CONTINUOUS, lb=0, name='beta_'+str(i))
    beta[0] = 0
    
    if args.quality:
        model.addConstrs((matches[t,i] - (alpha[t] if t in alpha else 0) - beta[i] - (beta[t-T] if t-T in beta else 0) <= 0 for t in range(1,T+K+1) if t in alpha or t-T in beta for i in beta  \
                if (t,i) in matches), "something...")
    else:
        model.addConstrs((COUNT((t,i)) - (alpha[t] if t in alpha else 0) - beta[i] - (beta[t-T] if t-T in beta else 0) <= 0 for t in range(1,T+K+1) if t in alpha or t-T in beta for i in beta \
                if (t,i) in matches), "something...")
    
    obj = quicksum(alpha[t] for t in alpha) + quicksum(beta[i] for i in beta)
    model.setObjective(obj, GRB.MINIMIZE)
    model.optimize()
<<<<<<< HEAD
    
    for t in range(T+1,T+K+1):
        if t not in alpha: continue
        if alpha[t].X > 0:
            print alpha[t]
    
    for i in range(1,K+1):
        results.append((demo[i+T-1], (beta[i].X if i in beta else 0)))
=======

    if args.graph_state:
        for i in range(1,K+1):
            results.append((demo[i+T-1], (beta[i].X if i in beta else 0), lpbeta[i]))
    else:
        for i in range(1,K+1):
            results.append((demo[i+T-1], (beta[i].X if i in beta else 0)))
>>>>>>> 6dfbeed68a5bb4446daaedd541ed3efc09fb4a55

if args.output:
    with open(args.output, 'w') as f:
        f.write(json.dumps(results))
else:
    print json.dumps(results)
