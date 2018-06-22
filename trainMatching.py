#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 20 11:21:22 2018

@author: kelseylieberman
"""
import argparse
import json
from gurobipy import *
import random
from sklearn import linear_model
from sklearn.preprocessing import PolynomialFeatures



parser = argparse.ArgumentParser()
parser.add_argument('--trainFiles', nargs = "+", help = "List of files to train on")
parser.add_argument('--testFiles', nargs = "+", help = "List of files to test")
parser.add_argument('-d', '--degree', default=1, type=int, help='type of polynomial to use while training')
parser.add_argument('-o', '--output')
parser.add_argument("-v", "--useVars", nargs = "+", type = int, default=[0, 1, 2, 3, 4, 5, 6, 7, 9, 11, 14, 15], \
                    help = "List of variables")
args = parser.parse_args()

results = ''
varsUsed = args.useVars
data = []
for fn in args.trainFiles:
    with open(fn, 'r') as f:
        data += json.load(f)
random.shuffle(data)

values = []
labels = []

for d in data:
    demo = d[0]
    values.append([demo[v] for v in varsUsed])
    labels.append(d[1])

poly = PolynomialFeatures(degree=args.degree)
X = poly.fit_transform(values)
LR = linear_model.LinearRegression()
LR.fit(X, labels)

for fn in args.testFiles:
    
    with open(fn, 'r') as f:
        data = json.load(f)

    K = data[0]
    T = data[1]
    matches = data[2]
    demo = data[3]
    
    testValues = [[demo[i][v] for v in varsUsed] for i in range(T,T+K)]
    X2 = poly.fit_transform(testValues)
    betaList = LR.predict(X2)
    beta = {i+1:betaList[i] for i in range(len(betaList))}
    beta[0] = 0
    
    quality = 0
    count = 0
    for t in range(T):
        values = {i:.1*random.random()+ matches[t][i] - beta[i] for i in beta if matches[t][i] != 0}
        max_i = max(values, key=values.get)
        if max_i != 0:
            count += 2
            quality += matches[t][max_i]
            del beta[max_i]
        else:
            count += 1
            quality += matches[t][max_i]
    
    
    model = Model('Online Matching')
    matchVars = {}
    for t in range(T,T+K):
        if t-T+1 not in beta: continue
        for i in beta:
            if matches[t][i] != 0:
                matchVars[(t,i)] = model.addVar(vtype = GRB.BINARY,  name = "match_" + str((t,i)))
    
    model.addConstrs((quicksum(matchVars[t,i] for i in range(1,K+1) if (t,i) in matchVars) <= 1 for t in range(T,T+K)), "only match with one other pair")
    model.addConstrs((quicksum(matchVars[t,i] for t in range(T,T+K) if (t,i) in matchVars) + quicksum(matchVars[i+T-1,j] for j in range(1, K+1) \
                               if (i+T-1,j) in matchVars) <= 1 for i in range(1,K+1)), "symmetry")
    
    obj = quicksum(matchVars[v]*matches[v[0]][v[1]] for v in matchVars)
    model.setObjective(obj, GRB.MAXIMIZE) 
    model.optimize()
    for v in matchVars:
        if matchVars[v].X != 0:
            count += 2
            quality += matches[v[0]][v[1]]
            
    results += str(count) + "\t" + str(quality) +"\n"
    
if args.output:
    with open(args.output, 'w') as f:
        f.write(results)
else:
    print results
    
    
    
    