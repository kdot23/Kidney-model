from gurobipy import *
import json
import pickle
import argparse
import os

parser = argparse.ArgumentParser(description="Computes the dual of our problem")
parser.add_argument('-i', '--input', default = 'data/', help='input directory to use')
parser.add_argument('--inputFile', help='input file to use, by default will use all files in directory')
parser.add_argument('-o', '--output', help='output file (json) to use. Includes demographic information and beta values')
parser.add_argument('--quality', action = "store_true", help="Optimize for quality")

args = parser.parse_args()

results = []
data = []

if args.inputFile:
    with open(args.input+'/'+args.inputFile, 'rb') as f:
        data.append(pickle.load(f))
else:
    for fn in os.listdir(args.input):
        with open(args.input+'/'+fn, 'rb') as f:
            data.append(pickle.load(f))

for d in data:
    #T is number of compatible pairs
    T = d[0]
    #K is number of incompatible pairs
    K = d[1]
    matches = d[2]
    demo = d[3]
    
    model = Model("Dual Optimizer")
    
    alpha = {}
    beta = {}
    
    for t in range(T+K):
        if any(k[0] == t for k in matches):
            alpha[t] = model.addVar(vtype = GRB.CONTINUOUS, lb=0, name='alpha_'+str(t))
    
    for i in range(1,K+1):
        if any(k[0] == i+T-1 for k in matches):
            beta[i] = model.addVar(vtype=GRB.CONTINUOUS, lb=0, name='beta_'+str(i))
        elif any(k[1] == i for k in matches):
            beta[i] = model.addVar(vtype=GRB.CONTINUOUS, lb=0, name='beta_'+str(i))
        elif any(k[2] == i for k in matches):
            beta[i] = model.addVar(vtype=GRB.CONTINUOUS, lb=0, name='beta_'+str(i))

    beta[0] = 0
    
    if args.quality:
        model.addConstrs((matches[t,i,j] - alpha[t] - beta[i] - beta[j] - (beta[t-T+1] if t+1-T in beta else 0) <= 0 for t in alpha for i in beta  \
                 for j in beta if (t,i,j) in matches), "something...")
    else:
        model.addConstrs((1 - alpha[t] - beta[i] - beta[j] - (beta[t-T+1] if t+1-T in beta else 0) <= 0 for t in alpha for i in beta \
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
