import json
import pickle
import argparse
import os
import pulp
from pulp import lpSum

parser = argparse.ArgumentParser(description="Computes the dual of our problem")
parser.add_argument('--inputFiles', nargs='+', default = ["data.dat"], help='input file to use')
parser.add_argument('-o', '--output', help='output file (json) to use. Includes demographic information and beta values')
parser.add_argument('--quality', action = "store_true", help="Optimize for quality")
parser.add_argument('--graph_state', action='store_true', help='Flag should be present if online LP estimation is included in training data')

args = parser.parse_args()

results = []

def COUNT(v):
    if v[1] == 0:
        return 1
    return 2

def calcBetasLP(C, matches, available_incompat):
    estimator = pulp.LpProblem("Estimate Betas", pulp.LpMinimize)
    posBetas = [i for i in available_incompat if any(k[0] == i+C for k in matches if k[1] in available_incompat) or \
            any(k[1] == i for k in matches if k[0] > C and k[0]-C in available_incompat)]
    beta = pulp.LpVariable.dicts('beta',posBetas,lowBound=0)
    estimator += lpSum(beta[i] for i in beta)
    if args.quality:
        for t in beta:
            for i in beta:
                if (t+C,i) not in matches: continue
                estimator += matches[t+C,i] - beta[i] -beta[t] <= 0, 'beta const ' + str((t,i))
    else:
        for t in beta:
            for i in beta:
                if (t+C,i) not in matches: continue
                estimator += COUNT((t+C,i)) - beta[i] -beta[t] <= 0, 'beta const ' + str((t,i))
    estimator.solve()
    newBeta =  {i:beta[i].value() for i in beta}
    for i in available_incompat:
        if i not in newBeta:
            newBeta[i] = 0
    return newBeta

for fn in args.inputFiles:
    with open(fn, 'rb') as f:
        d = pickle.load(f)
    #T is number of compatible pairs
    C = d[1]
    #K is number of incompatible pairs
    K = d[0]
    I = K
    T = d[2]
    matches = d[3]
    demo = d[5]
    departure_times = d[8]
    
    model = pulp.LpProblem('Calc Betas', pulp.LpMinimize)

    posAlphas = [t for t in range(1,C+1) if any(k[0] == t for k in matches)]
    posBetas = [i for i in range(1,K+1) if any(k[0] == i+C for k in matches) or any(k[1] == i for k in matches)]
    
    alpha = pulp.LpVariable.dicts('alpha',posAlphas,lowBound=0)
    beta = pulp.LpVariable.dicts('beta',posBetas,lowBound = 0)
    
    beta[0] = 0

    model += lpSum(alpha[t] for t in alpha) + lpSum(beta[i] for i in beta)
    if args.quality:
        for t in range(1,C+K+1):
            for i in beta:
                if (t,i) not in matches: continue
                model += matches[t,i] - (alpha[t] if t in alpha else 0) - beta[i] - (beta[t-C] if t-C in beta else 0) <= 0, 'alpha_beta '+str((t,i))
    else:
        for t in range(1,C+K+1):
            for i in beta:
                model += COUNT((t,i)) - (alpha[t] if t in alpha else 0) - beta[i] - (beta[t-C] if t-C in beta else 0) <= 0, 'alpha_beta '+str((t,i))

    
    
    model.solve()
    
    if args.graph_state:
        available_incompat = set()
        arriving_compat = {}
        arriving_incompat = {}
        for i in range(C):
            if demo[i][20] not in arriving_compat:
                arriving_compat[demo[i][20]] = []
            arriving_compat[demo[i][20]].append(i+1)
        for i in range(C,C+I):
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
            lpBeta = calcBetasLP(C, matches, available_incompat)
            for i in lpBeta:
                results.append((demo[i+C-1], (beta[i].value() if i in beta else 0), lpBeta[i]))
    else:
        for i in range(1,K+1):
            results.append((demo[i+C-1], (beta[i].value() if i in beta else 0)))

if args.output:
    with open(args.output, 'w') as f:
        f.write(json.dumps(results))
else:
    print json.dumps(results)
