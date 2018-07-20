#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 20 11:21:22 2018

Finds the optimal pair solutions using files that were trained by the online method.
Output is csv format of optimal count and quality for the population.
"""
import argparse
import json
import pickle
import random
from sklearn import linear_model
from sklearn import ensemble
from sklearn.preprocessing import PolynomialFeatures
import os
import pulp
from pulp import lpSum



parser = argparse.ArgumentParser()
parser.add_argument('--trainFiles', nargs = "+", help = "List of files to train on")
parser.add_argument('--testFiles', nargs = "+", help = "List of files to test")
parser.add_argument('--agents', help='output the quality of each agent to this file (.csv)')
parser.add_argument('-d', '--degree', default=1, type=int, help='type of polynomial to use while training')
parser.add_argument('-o', '--output', help = 'csv file to output count and quality to')
parser.add_argument("-v", "--useVars", nargs = "+", type = int, default=[0, 1, 2, 3, 4, 5, 6, 7, 9, 11, 14, 15], \
                    help = "List of variables")
parser.add_argument('--quality', action='store_true', help='Flag should be present if optimization should be done for quality')
parser.add_argument('--forestRegression', nargs='?', const=10, type=int, \
        help='Flag should be present if forest regression is to be used instead of Linear, optional argument of number of trees')
parser.add_argument('--graph', help='stem of output file for graph')
parser.add_argument('--lpEstimator', action='store_true', help='flag should be present if dual on incompatible pool only should be \
        used to estimate beta values')
parser.add_argument('--lpRepeat', action='store_true', help='flag should be present if lp repeat method is used to estimate betas')
parser.add_argument('-q', '--cadence',  type=int)
parser.add_argument('--incompatible_online', type=float, help='threshhold for probability an incompatible stays before it is matched')
parser.add_argument('--gamma', default=.9, type=float, help='gamma value used for calculating survival')
parser.add_argument('--graph_state', action='store_true', help='Flag should be present if lp estimator beta values is used for training')
args = parser.parse_args()


graph_colors = ["red", "blue", "green", "black"]

def COUNT(v):
    if v[1] == 0:
        return 1
    return 2

def getBloodTypes(demo):
    bd,br=0,0
    if demo[0]:
        br = 0
    elif demo[1]:
        br = 1
    elif demo[2]:
        br = 2
    elif demo[3]:
        br = 3
    if demo[4]:
        bd = 0
    elif demo[5]:
        bd = 1
    elif demo[6]:
        bd = 2
    elif demo[7]:
        bd = 3
    return br,bd

def calcBetasLP(C, matches, available_incompat):
    posBetas = [i for i in available_incompat if any(k[0] == i+C for k in matches if k[1] in available_incompat) or \
            any(k[1] == i for k in matches if k[0] > C and k[0]-C in available_incompat)]
    beta = pulp.LpVariable.dicts('beta', posBetas, lowBound = 0)
    posBetas = set(posBetas)
    estimator = pulp.LpProblem('Estimate Betas', pulp.LpMinimize)
    estimator += lpSum(beta[i] for i in posBetas)
    if args.quality:
        for t in posBetas:
            for i in posBetas:
                if (t+C,i) not in matches: continue
                estimator += matches[t+C,i] - beta[i] - beta[t] <= 0, 'beta qual ' + str((t,i))
    else:
        for t in posBetas:
            for i in posBetas:
                if (t+C,i) not in matches: continue
                estimator += COUNT((t+C,i)) - beta[i] - beta[t] <= 0, 'beta qual ' + str((t,i))
    estimator.solve()
    newBeta =  {i:beta[i].value() for i in posBetas}
    for i in available_incompat:
        if i not in newBeta:
            newBeta[i] = 0
    return newBeta


results = ''
graph = "digraph G {\n"
varsUsed = args.useVars
data = []
if not (args.lpEstimator or args.lpRepeat):
    for fn in args.trainFiles:
        with open(fn, 'r') as f:
            data += json.load(f)
    random.shuffle(data)

    values = []
    labels = []

    for d in data:
        demo = d[0]
        if args.graph_state:
            values.append([demo[v] for v in varsUsed] + [d[2]])
        else: 
            values.append([demo[v] for v in varsUsed])
        labels.append(d[1])

    poly = PolynomialFeatures(degree=args.degree)
    X = poly.fit_transform(values)
    if not args.forestRegression:
        LR = linear_model.Ridge()
    else:
        LR = ensemble.RandomForestRegressor(n_estimators=args.forestRegression)
    LR.fit(X, labels)

dataIndex = 0
agentInfo = ''
for fn in args.testFiles:
    
    with open(fn, 'r') as f:
        data = pickle.load(f)

    I = data[0]
    C = data[1]
    T = data[2]
    matches = data[3]
    demo = data[5]
    directed_matches = data[7]
    departure_times = data[8]
    available_incompat = set()
    unmatched_incompat = set()
    arriving_compat = {}
    arriving_incompat = {}
    num_rounds_present = {}
    lastBeta = {}
    count = 0
    quality = 0
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
                del num_rounds_present[i]
        available_incompat = available_incompat.difference(departing_incompat)
        unmatched_incompat = unmatched_incompat.union(set((i,lastBeta[i]) for i in departing_incompat))
        for i in available_incompat:
            num_rounds_present[i] += 1
        if t in arriving_incompat:
            available_incompat = available_incompat.union(arriving_incompat[t])
            for i in arriving_incompat[t]:
                num_rounds_present[i] = 0


        if args.lpEstimator or args.lpRepeat:
             #update betas one way
             beta = calcBetasLP(C, matches, available_incompat)
        else:
            #update betas another way
            beta = {}
            if args.graph_state:
                lpBeta = calcBetasLP(C, matches, available_incompat)
            for i in available_incompat:
                if args.graph_state:
                    testValue = [[demo[i+C-1][v] for v in varsUsed] + [lpBeta[i]]]
                else:
                    testValue = [[demo[i+C-1][v] for v in varsUsed]]
                testValue = poly.fit_transform(testValue)
                beta[i] = LR.predict(testValue)[0]
        beta[0] = 0
        lastBeta = beta

        if t in arriving_compat:
            #update betas

            #Do compatible mathcing stuff
            for i in arriving_compat[t]:
                if args.quality:
                    values = {j:matches[i,j]-beta[j] for j in available_incompat.union(set([0])) if (i,j) in matches }
                else:
                    values = {j:COUNT((i,j))-beta[j] for j in available_incompat.union(set([0])) if (i,j) in matches}
                max_index = max(values, key=values.get)
                if max_index == 0:
                    count += 1
                    agentInfo += "C" + str(i) + "\t" + str(t) + "\t" + str(directed_matches[i,0]) + "\t" \
                    + "C" + "\t" + str(directed_matches[i,0]) + "\t" + "C" + "\t" + str(demo[i-1][20]) + "\t" + str(t) + "\t" + str(0) + "\n"
                else:
                    available_incompat.remove(max_index)
                    count += 2
                    agentInfo += "C" + str(i) + "\t" + str(t) + "\t" + str(directed_matches[max_index+C,i]) + "\t" \
                    + "I" + "\t" + str(directed_matches[i,max_index+C]) + "\t" + "I" + str(demo[i-1][20]) + "\t" + str(t) + "\t" + str(0) + "\n"
                    agentInfo += "I" + str(max_index) + "\t" + str(t) + "\t" + str(directed_matches[i,max_index+C]) + "\t" \
                    + "C" + "\t" + str(directed_matches[max_index+C,i]) + "\t" + "C" + "\t" + str(demo[max_index+C-1][20]) + "\t" + \
                    str(departure_times[max_index-1]) + "\t" + str(str(beta[max_index]) + "\n"
                    if args.lpRepeat:
                        beta = calcBetasLP(C, matches, available_incompat)
                        beta[0] = 0
                quality += matches[i,max_index]


        if args.incompatible_online:
            probs = {i:demo[i+C-1][19]*args.gamma**num_rounds_present[i] for i in available_incompat}
            for i in probs:
                if i not in available_incompat: continue
                if probs[i] < args.incompatible_online:
                    values = {j:matches[i+C,j] - beta[j] for j in available_incompat if (i+C,j) in matches}
                    if len(values) == 0: continue
                    max_index = max(values, key=values.get)
                    if values[max_index] > 0:
                        available_incompat.remove(i)
                        available_incompat.remove(max_index)
                        quality += matches[i+C,max_index]
                        count += COUNT((i+C,max_index))
                        agentInfo += "I" + str(i) + "\t" + str(t) + "\t" + str(directed_matches[max_index+C,i+C]) + "\t" \
                        + "I" + "\t" + str(directed_matches[i+C,max_index+C]) + "\t" + "I" + "\t" + \
                        str(demo[i-1][20]) + "\t" + str(departure_times[i-C-1]) + "\t" + str(beta[i]) + "\n"
                        agentInfo += "I" + str(max_index) + "\t" + str(t) + "\t" + str(directed_matches[i+C,max_index+C]) + "\t" \
                        + "I" + "\t" + str(directed_matches[max_index+C,i+C]) + "\t" + "I" + "\t" +  \
                        str(demo[max_index+C-1][20]) + "\t" + str(departure_times[max_index-1]) + "\t" + str(beta[max_index]) + "\n"
        if args.cadence and (t)%args.cadence==0:
            #Do incompatible matching stuff
            model = pulp.LpProblem('incompatible matching', pulp.LpMaximize)
            matchVars = [(i+C,j) for i in available_incompat for j in available_incompat if (i+C,j) in matches]
            x = pulp.LpVariable.dicts('match', matchVars, lowBound = 0, upBound = 1, cat = pulp.LpInteger)
            matchVars = set(matchVars)
            if args.quality:
                model += lpSum(x[v]*matches[v] for v in matchVars)
            else:
                model += lpSum(x[v]*COUNT(v) for v in matchVars)
            for t in available_incompat:
                model += lpSum(x[t+C,i] for i in available_incompat if (t+C,i) in matchVars) <= 1, 'only match with 1 '+str(t)
            for i in available_incompat:
                model += lpSum(x[t+C,i] for t in available_incompat if (t+C,i) in matchVars) + lpSum(x[i+C,j] for j in available_incompat \
                        if (i+C,j) in matchVars) <= 1, 'symetry '+str(i)
            model.solve()
            count += sum(COUNT(v)*x[v].value() for v in matchVars)
            quality += sum(matches[v]*x[v].value() for v in matchVars)

            for v in matchVars:
                if round(x[v].value()) != 0:
                    available_incompat.remove(v[0]-C)
                    available_incompat.remove(v[1])
                    #Agent Info Stuff
                    agentInfo += "I" + str(v[0]-C) + "\t" + str(t) + "\t" + str(directed_matches[v[1]+C,v[0]]) + "\t" \
                    + "I" + "\t" + str(directed_matches[v[0],v[1]+C]) + "\t" + "I" + "\t" + \
                    str(demo[v[0]-1][20]) + "\t" + str(departure_times[v[0]-C-1]) + "\t" + str(beta[v[0]-C]) + "\n"
                    agentInfo += "I" + str(v[1]) + "\t" + str(t) + "\t" + str(directed_matches[v[0],v[1]+C]) + "\t" \
                    + "I" + "\t" + str(directed_matches[v[1]+C,v[0]]) + "\t" + "I" + "\t" \
                    + str(demo[v[1]+C-1][20]) + "\t" + str(departure_times[v[1]-1]) + "\t" + str(beta[v[1]]) + "\n"
    unmatched_incompat = unmatched_incompat.union(set((i,lastBeta[i]) for i in available_incompat))
    for a in unmatched_incompat:
        i = a[0]
        b = a[1]
        agentInfo += "I" + str(i) + "\t" + str(T) + "\t" + str(0) + "\t" \
        + "N" + "\t" + str(0) + "\t" + "N" + "\t" + \
        + str(demo[i+C-1][20]) + "\t" + str(departure_times[i-1]) + "\t" + str(b) + "\n"

    results += str(count) + '\t' + str(quality) + '\n'




    
    """
    testValues = [[demo[i][v] for v in varsUsed] for i in range(T,T+K)]
    X2 = poly.fit_transform(testValues)
    if not (args.lpEstimator or args.lpRepeat):
        betaList = LR.predict(X2)
        beta = {i+1:betaList[i] for i in range(len(betaList))}
    else:
        beta = calcBetasLP(T,K,matches,used_incompat)
    
    beta[0] = 0
    
    quality = 0
    count = 0
    for t in range(1,T+1):
        if args.quality: 
            values = {i:.1*random.random() + matches[t, i] - beta[i] for i in beta if (t,i) in matches }
        else:
            values = {i:.1*random.random()+COUNT((t,i)) - beta[i] for i in beta if (t,i) in matches}
        max_index = max(values, key=values.get)
        
        if max_index == 0:
            agentInfo += "C" + str(t) + "\t" + str(t) + "\t" + str(directed_matches[t,0]) + "\t" \
           + "C" + "\t" + str(directed_matches[t,0]) + "\t" + "C" + "\t" + str(0) + "\n"
        else:
            agentInfo += "C" + str(t) + "\t" + str(t) + "\t" + str(directed_matches[max_index+T,t]) + "\t" \
           + "I" + "\t" + str(directed_matches[t,max_index+T]) + "\t" + "I" + "\t" + str(0) + "\n"
            agentInfo += "I" + str(max_index) + "\t" + str(t) + "\t" + str(directed_matches[t,max_index+T]) + "\t" \
           + "C" + "\t" + str(directed_matches[max_index+T,t]) + "\t" + "C" + "\t" + str(beta[max_index]) + "\n"
        if max_index != 0:
            count += 2
            quality += matches[t,max_index]
            used_incompat.add(max_index)
            del beta[max_index]
            if args.lpRepeat:
                beta = calcBetasLP(T, K, matches, used_incompat)
                beta[0] = 0
            bt1 = getBloodTypes(demo[t-1])
            bt2 = getBloodTypes(demo[max_index + T - 1])
            graph += "edge [color="+graph_colors[bt1[1]] + "];\n"
            graph += "C" + str(t) + " [color="+graph_colors[bt1[0]]+"];\n"
            graph += "I" + str(max_index) + " [color="+graph_colors[bt2[0]]+"];\n"
            graph += "C" + str(t) + " -> I" + str(max_index) + ";\n"
            graph += "edge [color="+graph_colors[bt2[1]] + "];\n"
            graph += "I" + str(max_index) + " -> C" + str(t) + ";\n"
        else:
            count += 1
            quality += matches[t,max_index]
            bt = getBloodTypes(demo[t-1])
            graph += "edge [color="+graph_colors[bt[1]] + "];\n"
            graph += "node [color="+graph_colors[bt[0]]+"];\n"
            graph += "C" + str(t) + " -> C" + str(t) + ";\n"
    
    
    model = Model('Online Matching')
    matchVars = {}
    for t in range(T+1,T+K+1):
        if t-T not in beta: continue
        for i in beta:
            if (t,i) in matches:
                matchVars[(t,i)] = model.addVar(vtype = GRB.BINARY,  name = "match_" + str((t,i)))
    
    model.addConstrs((quicksum(matchVars[t,i] for i in range(1,K+1) if (t,i) in matchVars) <= 1 for t in range(T,T+K+1)), "only match with one other pair")
    model.addConstrs((quicksum(matchVars[t,i] for t in range(T,T+K+1) if (t,i) in matchVars) + quicksum(matchVars[i+T,j] for j in range(1, K+1) \
                               if (i+T,j) in matchVars) <= 1 for i in range(1,K+1)), "symmetry")
    
    if args.quality:
        obj = quicksum(matchVars[v]*matches[v] for v in matchVars)
    else:
        obj = quicksum(matchVars[v]*COUNT(v) for v in matchVars)
    model.setObjective(obj, GRB.MAXIMIZE) 
    model.optimize()
    for v in matchVars:
        if round(matchVars[v].X) != 0:
            agentInfo += "I" + str(v[0]-T) + "\t" + str(T+1) + "\t" + str(directed_matches[v[1]+T,v[0]]) + "\t" \
            + "I" + "\t" + str(directed_matches[v[0],v[1]+T]) + "\t" + "I" + "\t" + str(beta[v[0]-T]) + "\n"
            agentInfo += "I" + str(v[1]) + "\t" + str(T+1) + "\t" + str(directed_matches[v[0],v[1]+T]) + "\t" \
            + "I" + "\t" + str(directed_matches[v[1]+T,v[0]]) + "\t" + "I" + "\t" + str(beta[v[1]]) + "\n"
     
            bt1 = getBloodTypes(demo[v[0]-1])
            bt2 = getBloodTypes(demo[v[1] + T - 1])
            graph += "edge [color="+graph_colors[bt1[1]] + "];\n"
            graph += "I" + str(v[0]-T) + " [color="+graph_colors[bt1[0]]+"];\n"
            graph += "I" + str(v[1]-1) + " [color="+graph_colors[bt2[0]]+"];\n"
            graph += "I" + str(v[0]-T) + " -> I" + str(v[1]-1) + ";\n"
            graph += "edge [color="+graph_colors[bt2[1]] + "];\n"
            graph += "node [color="+graph_colors[bt2[0]]+"];\n"
            graph += "I" + str(v[1]-1) + " -> I" + str(v[0]-T) + ";\n"
            del beta[v[0]-T]
            del beta[v[1]]
    del beta[0]       
    for i in beta:
        agentInfo += "I" + str(i) + "\t" + str(T+2) + "\t" + str(0) + "\t" \
        + "N" + "\t" + str(0) + "\t" + "N" + "\t" + str(beta[i]) + "\n"
    
    quality += sum(matchVars[v].X*matches[v] for v in matchVars)
    count += sum(COUNT(v)*matchVars[v].X for v in matchVars)

    results += str(count) + "\t" + str(quality) +"\n"
    graph += "}"
    if args.graph:
        with open(args.graph+str(dataIndex)+'.gv', 'w') as f:
            f.write(graph)
        os.system('dot -Tpdf ' + args.graph + str(dataIndex) +'.gv -o ' + args.graph + str(dataIndex) + '.pdf')
    dataIndex += 1
    """
    
if args.output:
    with open(args.output, 'w') as f:
        f.write(results)
else:
    print results
    
if args.agents:
    with open(args.agents, 'w') as f:
        f.write(agentInfo)   
