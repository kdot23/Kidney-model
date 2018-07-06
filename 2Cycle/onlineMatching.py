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
from gurobipy import *
import random
from sklearn import linear_model
from sklearn import ensemble
from sklearn.preprocessing import PolynomialFeatures
import os



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

results = ''
graph = "digraph G {\n"
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
if not args.forestRegression:
    LR = linear_model.LinearRegression()
else:
    LR = ensemble.RandomForestRegressor(n_estimators=args.forestRegression)
LR.fit(X, labels)

dataIndex = 0
agentInfo = ''
for fn in args.testFiles:
    
    with open(fn, 'r') as f:
        data = pickle.load(f)

    K = data[0]
    T = data[1]
    matches = data[2]
    demo = data[4]
    directed_matches = data[6]
    
    testValues = [[demo[i][v] for v in varsUsed] for i in range(T,T+K)]
    X2 = poly.fit_transform(testValues)
    betaList = LR.predict(X2)
    beta = {i+1:betaList[i] for i in range(len(betaList))}
    
    beta[0] = 0
    
    quality = 0
    count = 0
    for t in range(1,T+1):
        if args.quality: 
            values = {i:.1*random.random()+ matches[t,i] - beta[i] for i in beta if (t,i) in matches }
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
            del beta[max_index]
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
    
    if (args.incompatibleOnly):
        quality = sum(matchVars[v].X*matches[v] for v in matchVars) + sum(matches[(v,0)] for v in range(T))
        count = sum(COUNT(v)*matchVars[v].X for v in matchVars) + T
    else:
        quality = sum(matchVars[v].X*matches[v] for v in matchVars)
        count = sum(COUNT(v)*matchVars[v].X for v in matchVars)

    results += str(count) + "\t" + str(quality) +"\n"
    graph += "}"
    if args.graph:
        with open(args.graph+str(dataIndex)+'.gv', 'w') as f:
            f.write(graph)
        os.system('dot -Tpdf ' + args.graph + str(dataIndex) +'.gv -o ' + args.graph + str(dataIndex) + '.pdf')
    dataIndex += 1
    
if args.output:
    with open(args.output, 'w') as f:
        f.write(results)
else:
    print results
    
if args.agents:
    with open(args.agents, 'w') as f:
        f.write(agentInfo)   
