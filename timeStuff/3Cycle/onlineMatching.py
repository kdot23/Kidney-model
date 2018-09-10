#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 20 11:21:22 2018

Finds the optimal pair solutions using files that were trained by the online method.
Output is csv format of optimal count and quality for the population.
"""
import argparse
import json
import random
from sklearn import linear_model
from sklearn import ensemble
from sklearn.preprocessing import PolynomialFeatures
import os
import pulp
from pulp import lpSum
from zipfile import ZipFile



parser = argparse.ArgumentParser()
parser.add_argument('--trainFiles', nargs = "+", help = "List of files to train on")
parser.add_argument('--trainZipFile', help='filename of zip file to look for trainFiles in, if not given data is assumed to be uncompressed')
parser.add_argument('--testFiles', nargs = "+", help = "List of files to test")
parser.add_argument('--testZipFile', help='filename of zip file to look for testFiles in, if not given data is assumed to be uncompressed')
parser.add_argument('-o', '--output', help = 'csv file to output count and quality to')
parser.add_argument('--agents', help='output the quality of each agent to this file (.csv)')
parser.add_argument("-v", "--useVars", nargs = "+", type = int, default=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18], \
                    help = "List of variables")
parser.add_argument('--quality', action='store_true', help='Flag should be present if optimization is being done for quality')
parser.add_argument('--forestRegression', nargs='?', const=10, type=int, \
        help='Flag should be present if forest regression is to be used instead of Linear, optional argument of number of trees')
parser.add_argument('--graph', help='stem of output file for graph')
parser.add_argument('-d', '--degree', default=1, type=int, help='type of polynomial to use while training')
parser.add_argument('--lpEstimator', action='store_true', help='should be present if dual on incompatibles only should be used to \
        estimate betas')
parser.add_argument('--lpRepeat', action='store_true', help='should be present if lp repeat method is used to estimate betas')
parser.add_argument('-q', '--cadence', default=1, type=int, help='frequency to clear incompatible pool')
parser.add_argument('--incompatible_online', action='store_true', help='threshhold for probability an incompatible stays before it is matched')
parser.add_argument('--gamma', default=.9, type=float, help='gamma value used for calculating survival')
parser.add_argument('--graph_state', action='store_true', help='Flag should be present if online LP estimation is included in training data')
args = parser.parse_args()


graph_colors = ["red", "blue", "green", "black"]


def convertStringToTuple(s):
    return tuple(int(i) for i in s.split(','))

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

def COUNT(v):
    if v[1] == 0:
        return 1
    if v[2] == 0:
        return 2
    return 3
def calcBetaLP(C, matches, available_incompat):
    estimator = pulp.LpProblem('estimate betas', pulp.LpMinimize)
    posBetas = [i for i in available_incompat if any(k[0] == i+T for k in matches if k[1] in available_incompat and k[2] in available_incompat) or \
            any(k[1] == i for k in matches if k[0] > C and k[0]-C in available_incompat and k[2] in available_incompat) or \
            any(k[2] == i for k in matches if k[0] > C and k[0]-C in available_incompat and k[1] in available_incompat)]
    beta = pulp.LpVariable.dicts('beta',posBetas,lowBound=0)
    posBetas = set(posBetas)
    estimator += lpSum(beta[i] for i in posBetas)
    for t in posBetas:
        for i in posBetas:
            for j in posBetas:
                if (t+C,i,j) not in matches: continue
                if args.quality:
                    estimator += matches[t+C,i,j] - beta[i] -beta[j] -beta[t] <= 0, 'beta constraint ' + str((t,i,j))
                else:
                    estimator += COUNT((t+C,i,j)) - beta[i] -beta[j] -beta[t] <= 0, 'beta constraint ' + str((t,i,j))
    estimator.solve()
    newBeta = {i:beta[i].value() for i in posBetas}
    for i in available_incompat:
        if i not in newBeta:
            newBeta[i] = 0
    return newBeta

results = ''
agentInfo = '' 
#Pair id, time of match, get quality, get donor type, give quality, give recipient type, beta value (0 for incompatible)
#time of match = T+1 for incompatible pairs matched at the end
graph = "digraph G {\n"
varsUsed = args.useVars
data = []
if args.trainFiles:
    if args.trainZipFile:
        trainZipFile = ZipFile(args.trainZipFile)
    for fn in args.trainFiles:
        if args.trainZipFile:
            data += json.loads(trainZipFile.read(fn))
        else:
            with open(fn, 'r') as f:
                data += json.load(f)
    if args.trainZipFile:
        trainZipFile.close()
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
if args.testZipFile:
    testZipFile = ZipFile(args.testZipFile)
for fn in args.testFiles:
    if args.testZipFile:
        data = json.loads(testZipFile.read(fn))
    else:
        with open(fn, 'r') as f:
            data = json.load(f)

    I = data[0]
    C = data[1]
    T = data[2]
    matches = data[4]
    matches = {convertStringToTuple(i):matches[i] for i in matches}
    demo = data[5]
    directed_matches = data[7]
    directed_matches = {convertStringToTuple(i):directed_matches[i] for i in directed_matches}
    departure_times = data[8]
    used_incompat = set()

    available_incompat = set()
    unmatched_incompat = set()
    arriving_incompat = {}
    arriving_compat = {}
    num_rounds_present = {}

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

    lastBeta = {}
    for t in range(T+1):
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
            beta = calcBetaLP(C, matches, available_incompat)
        else:
            beta = {}
            if args.graph_state:
                lpBeta = calcBetaLP(C, matches, available_incompat)
            for i in available_incompat:
                if args.graph_state:
                    testValue = [[demo[i-1][v] for v in varsUsed] + [lpBeta[i]]]
                else:
                    testValue = [[demo[i-1][v] for v in varsUsed]]
                testValue = poly.fit_transform(testValue)
                beta[i] = LR.predict(testValue)[0]
        beta[0] = 0
        lastBeta = beta

        if t in arriving_compat:
            for i in arriving_compat[t]:
                if args.quality:
                    values = {(i,j,k):matches[i,j,k] - beta[j] - beta[k] for j in available_incompat.union(set([0]))
                            for k in available_incompat.union(set([0])) if (i,j,k) in matches}
                else:
                    values = {(i,j,k):COUNT((i,j,k)) - beta[j] - beta[k] for j in available_incompat.union(set([0]))
                            for k in available_incompat.union(set([0])) if (i,j,k) in matches}
                max_index = max(values,key=values.get)
                quality += matches[max_index]
                if max_index[1] == 0:
                    count += 1
                    agentInfo += "C" + str(i) + "\t" + str(t) + "\t" + str(directed_matches[max_index[0],0]) + "\t" \
                    + "C" + "\t" + str(directed_matches[max_index[0],0]) + "\t" + "C" + "\t"   + str(demo[max_index[0]-1][20]) + "\t" + str(t)+ "\t" + str(0) + "\n"
                elif max_index[2] == 0:
                    count += 2
                    available_incompat.remove(max_index[1])
                    agentInfo += "C" + str(i) + "\t" + str(t) + "\t" + str(directed_matches[max_index[1]+C,max_index[0]]) + "\t" \
                    + "I" + "\t" + str(directed_matches[max_index[0],max_index[1]+C]) + "\t" + "I" + "\t"   + str(demo[max_index[0]-1][20]) + "\t" + str(t) + "\t" + str(0) +"\n"
                    agentInfo += "I" + str(max_index[1]) + "\t" + str(t) + "\t" + str(directed_matches[max_index[0],max_index[1]+C]) + "\t" \
                    + "C" + "\t" + str(directed_matches[max_index[1]+C,max_index[0]]) + "\t" + "C" + "\t"   + str(demo[max_index[1]+C-1][20]) + "\t" \
                    + str(departure_times[max_index[1]-1])+ "\t" + str(beta[max_index[1]]) + "\n"
                else:
                    count += 3
                    available_incompat.remove(max_index[1])
                    available_incompat.remove(max_index[2])
                    agentInfo += "C" + str(i) + "\t" + str(t) + "\t" + str(directed_matches[max_index[2]+C,max_index[0]]) + "\t" \
                    + "I" + "\t" + str(directed_matches[max_index[0],max_index[1]+C]) + "\t" + "I"   +  "\t" + str(demo[max_index[0]-1][20]) + "\t" + str(t)+ "\t"+ str(0) +"\n"
                    agentInfo += "I" + str(max_index[1]) + "\t" + str(t) + "\t" + str(directed_matches[max_index[0],max_index[1]+C]) + "\t" \
                    + "C" + "\t" + str(directed_matches[max_index[1]+C,max_index[2]+C]) + "\t" + "I" + "\t"   + str(demo[max_index[1]+C-1][20]) \
                    + "\t" + str(departure_times[max_index[1]-1])+ "\t"+ str(beta[max_index[1]]) +"\n"
                    agentInfo += "I" + str(max_index[2]) + "\t" + str(t) + "\t" + str(directed_matches[max_index[1]+C,max_index[2]+C]) + "\t" \
                    + "I" + "\t" + str(directed_matches[max_index[2]+C,max_index[0]]) + "\t" + "C"   + "\t" + str(demo[max_index[2]+C-1][20]) \
                    + "\t" + str(departure_times[max_index[2]-1])+ "\t"+ str(beta[max_index[2]]) +"\n"

                if args.lpRepeat:
                    beta = calcBetaLP(C, matches, available_incompat)
                    beta[0] = 0


        if args.incompatible_online:
            matching_incompat = set(i for i in available_incompat if departure_times[i-1] <= t)
            for i in matching_incompat:
                if i not in available_incompat: continue
                values = {(i+C,j,k):matches[i+C,j, k] - beta[j] - beta[k] for j in available_incompat for k in available_incompat.union(set([0])) if (i+C,j,k) in matches}
                if len(values) == 0: continue
                max_index = max(values, key=values.get)
                available_incompat.remove(i)
                available_incompat.remove(max_index[1])
                quality += matches[max_index]
                count += COUNT(max_index)
                if max_index[2] == 0:
                    agentInfo += "I" + str(i) + "\t" + str(t) + "\t" + str(directed_matches[max_index[1]+C,max_index[0]]) + "\t" \
                    + "I" + "\t" + str(directed_matches[max_index[0],max_index[1]+C]) + "\t" + "I" + "\t" + str(demo[i+C-1][20]) +\
                    '\t' + str(departure_times[i-1]) + '\t' +  str(beta[max_index[0]-C]) + "\n"
                    agentInfo += "I" + str(max_index[1]) + "\t" + str(t) + "\t" + str(directed_matches[max_index[0],max_index[1]+C]) + "\t" \
                    + "I" + "\t" + str(directed_matches[max_index[1]+C,max_index[0]]) + "\t" + "I" + "\t" + str(demo[max_index[1]+C-1][20]) \
                    + '\t' + str(departure_times[max_index[1]-1]) + '\t' + str(beta[max_index[1]]) + "\n"
                else:
                    available_incompat.remove(max_index[2])
                    agentInfo += "I" + str(i) + "\t" + str(t) + "\t" + str(directed_matches[max_index[2]+C,max_index[0]]) + "\t" \
                    + "I" + "\t" + str(directed_matches[max_index[0],max_index[1]+C]) + "\t" + "I" + "\t" +\
                    str(demo[i+C-1][20]) + '\t' + str(departure_times[i-1]) +'\t' + str(beta[max_index[0]-C]) + "\n"
                    agentInfo += "I" + str(max_index[1]) + "\t" + str(t) + "\t" + str(directed_matches[max_index[0],max_index[1]+C]) + "\t" \
                    + "I" + "\t" + str(directed_matches[max_index[1]+C,max_index[2]+C]) + "\t" + "I" + "\t" +\
                    str(demo[max_index[1]+C-1][20]) + '\t' + str(departure_times[max_index[1]-1]) + '\t' +  str(beta[max_index[1]]) + "\n"
                    agentInfo += "I" + str(max_index[2]) + "\t" + str(t) + "\t" + str(directed_matches[max_index[1]+C,max_index[2]+C]) + "\t" \
                    + "I" + "\t" + str(directed_matches[max_index[2]+C,max_index[0]]) + "\t" + "I" + "\t" + \
                    str(demo[max_index[2]+C-1][20]) + '\t' + str(departure_times[max_index[2]-1]) + '\t' + str(beta[max_index[2]]) + "\n"
        if args.cadence and t%args.cadence == 0:
            model = pulp.LpProblem('match incompatibles', pulp.LpMaximize)
            matchVars = [(i+C,j,k) for i in available_incompat for j in available_incompat for k in available_incompat if (i+C,j,k) in matches]
            x = pulp.LpVariable.dicts('match',matchVars,lowBound=0,upBound=1,cat=pulp.LpInteger)
            matchVars = set(matchVars)
            if args.quality:
                model += lpSum(x[v]*matches[v] for v in matchVars)
            else:
                model += lpSum(x[v]*COUNT(v) for v in matchVars)
            for t in available_incompat:
                model += lpSum(x[t+C,i,j] for i in available_incompat for j in available_incompat if (t+C,i,j) in matchVars) <= 1,\
                        'only match with one '+str(t)
            for i in available_incompat:
                model += lpSum(x[t+C,i,j] for t in available_incompat for j in available_incompat if (t+C,i,j) in matchVars) + \
                        lpSum(x[t+C,j,i] for t in available_incompat for j in available_incompat if (t+C,j,i) in matchVars) + \
                        lpSum(x[i+C,j,k] for j in available_incompat for k in available_incompat if (i+C,j,k) in matchVars) <= 1, \
                        'symetry '+str(i)
            model.solve()
            count += sum(COUNT(v)*x[v].value() for v in matchVars)
            quality += sum(matches[v]*x[v].value() for v in matchVars)

            for v in matchVars:
                if round(x[v].value()) != 0:
                    available_incompat.remove(v[0]-C)
                    available_incompat.remove(v[1])
                    if v[2] == 0:
                        agentInfo += "I" + str(v[0]-C) + "\t" + str(t) + "\t" + str(directed_matches[v[1]+C,v[0]]) + "\t" \
                        + "I" + "\t" + str(directed_matches[v[0],v[1]+C]) + "\t" + "I"   + "\t" + str(demo[v[0]-1][20]) + "\t" + str(departure_times[v[0]-C-1])+ "\t"+ str(beta[v[0]-C]) +"\n"
                        agentInfo += "I" + str(v[1]) + "\t" + str(t) + "\t" + str(directed_matches[v[0],v[1]+C]) + "\t" \
                        + "I" + "\t" + str(directed_matches[v[1]+C,v[0]]) + "\t" + "I" + "\t"   + str(demo[v[1]+C-1][20]) + "\t" + str(departure_times[v[1]-1])+ "\t"+ str(beta[v[1]]) +"\n"
                    else:
                        available_incompat.remove(v[2])
                        agentInfo += "I" + str(v[0]-C) + "\t" + str(t) + "\t" + str(directed_matches[v[2]+C,v[0]]) + "\t" \
                        + "I" + "\t" + str(directed_matches[v[0],v[1]+C]) + "\t" + "I" + "\t"   + str(demo[v[0]-1][20]) + "\t" + str(departure_times[v[0]-C-1])+ "\t"+ str(beta[v[0]-C]) + "\n"
                        agentInfo += "I" + str(v[1]) + "\t" + str(t) + "\t" + str(directed_matches[v[0],v[1]+C]) + "\t" \
                        + "I" + "\t" + str(directed_matches[v[1]+C,v[2]+C]) + "\t" + "I" + "\t"   \
                        + str(demo[v[1]+C-1][20]) + "\t" + str(departure_times[v[1]-1])+ "\t"+ str(beta[v[1]]) +"\n"
                        agentInfo += "I" + str(v[2]) + "\t" + str(t) + "\t" + str(directed_matches[v[1]+C,v[2]+C]) + "\t" \
                                + "I" + "\t" + str(directed_matches[v[2]+C,v[0]]) + "\t" + "I"  + "\t" \
                                + str(demo[v[2]+C-1][20]) + "\t" + str(departure_times[v[2]-1])+ "\t" + str(beta[v[2]]) +"\n"

    results += str(count) + '\t' + str(quality) + '\n'
    unmatched_incompat = unmatched_incompat.union(set((i,lastBeta[i]) for i in available_incompat))
    for a in unmatched_incompat:
        i = a[0]
        b = a[1]
        agentInfo += "I" + str(i) + "\t" + str(T) + "\t" + str(0) + "\t" + "N" + "\t" + str(0) + "\t" + "N"  +"\t" + str(demo[i+C-1][20]) + "\t" + str(departure_times[i-1])+ "\t" + str(b) + "\n"

if args.testZipFile:
    testZipFile.close()
    """


    
    testValues = [[demo[i][v] for v in varsUsed] for i in range(T,T+K)]
    X2 = poly.fit_transform(testValues)
    if not (args.lpEstimator or args.lpRepeat):
        betaList = LR.predict(X2)
        beta = {i+1:betaList[i] for i in range(len(betaList))}
    else:
        beta = calcBetaLP(T, K, matches, used_incompat)

    #FOR CAPPING
    for i in beta:
        if beta[i] < 0:
            beta[i] = 0  
    beta[0] = 0
    
    quality = 0
    count = 0
    for t in range(1,T+1):
        values = {(t, i, j):.1*random.random()+ matches[t,i,j] - beta[i] - beta[j] \
                  for i in beta for j in beta if (t,i,j) in matches}
        max_index = max(values, key=values.get)
        count += COUNT(max_index)
        quality += matches[ax_index]
       
        if max_index[1] == 0:
            agentInfo += "C" + str(t) + "\t" + str(t) + "\t" + str(directed_matches[max_index[0],0]) + "\t" \
           + "C" + "\t" + str(directed_matches[max_index[0],0]) + "\t" + "C" + "\t" + str(0) + "\n"
        elif max_index[2] == 0:
            agentInfo += "C" + str(t) + "\t" + str(t) + "\t" + str(directed_matches[max_index[1]+T,max_index[0]]) + "\t" \
           + "I" + "\t" + str(directed_matches[max_index[0],max_index[1]+T]) + "\t" + "I" + "\t" + str(0) + "\n"
            agentInfo += "I" + str(max_index[1]) + "\t" + str(t) + "\t" + str(directed_matches[max_index[0],max_index[1]+T]) + "\t" \
           + "C" + "\t" + str(directed_matches[max_index[1]+T,max_index[0]]) + "\t" + "C" + "\t" + str(beta[max_index[1]]) + "\n"
        else:
            agentInfo += "C" + str(t) + "\t" + str(t) + "\t" + str(directed_matches[max_index[2]+T,max_index[0]]) + "\t" \
           + "I" + "\t" + str(directed_matches[max_index[0],max_index[1]+T]) + "\t" + "I" + "\t" + str(0) +  "\n"
            agentInfo += "I" + str(max_index[1]) + "\t" + str(t) + "\t" + str(directed_matches[max_index[0],max_index[1]+T]) + "\t" \
           + "C" + "\t" + str(directed_matches[max_index[1]+T,max_index[2]+T]) + "\t" + "I" + "\t" + str(beta[max_index[1]]) + "\n"
            agentInfo += "I" + str(max_index[2]) + "\t" + str(t) + "\t" + str(directed_matches[max_index[1]+T,max_index[2]+T]) + "\t" \
           + "I" + "\t" + str(directed_matches[max_index[2]+T,max_index[0]]) + "\t" + "C" + "\t" + str(beta[max_index[2]]) + "\n"

        if max_index[1] != 0:
            del beta[max_index[1]]
            used_incompat.add(max_index[1])
        if max_index[2] != 0:
            del beta[max_index[2]] 
            used_incompat.add(max_index[2])
        if args.lpRepeat and (max_index[2] != 0 or max_index[1] != 0):
            beta = calcBetaLP(T, K, matches, used_incompat)
            beta[0] = 0
            

    """
    """
        if max_i != 0:
            count += 2
            quality += matches[t][max_i]
            del beta[max_i]
            bt1 = getBloodTypes(demo[t])
            bt2 = getBloodTypes(demo[max_i + T - 1])
            graph += "edge [color="+graph_colors[bt1[1]] + "];\n"
            graph += "C" + str(t) + " [color="+graph_colors[bt1[0]]+"];\n"
            graph += "I" + str(max_i-1) + " [color="+graph_colors[bt2[0]]+"];\n"
            graph += "C" + str(t) + " -> I" + str(max_i-1) + ";\n"
            graph += "edge [color="+graph_colors[bt2[1]] + "];\n"
            graph += "I" + str(max_i-1) + " -> C" + str(t) + ";\n"
        else:
            count += 1
            quality += matches[t][max_i]
            bt = getBloodTypes(demo[t])
            graph += "edge [color="+graph_colors[bt[1]] + "];\n"
            graph += "node [color="+graph_colors[bt[0]]+"];\n"
            graph += "C" + str(t) + " -> C" + str(t) + ";\n"
    """
    """
         
             
    model = Model('Online Matching')
    matchVars = {}
    for t in range(T+1,T+K+1):
        if t-T not in beta: continue
        for i in beta:
            for j in beta:
                if (t,i,j) in matches:
                    matchVars[(t,i,j)] = model.addVar(vtype = GRB.BINARY,  name = "match_" + str((t,i,j)))
    
    model.addConstrs((quicksum(matchVars[t,i,j] for i in range(K+1) for j in range(K+1) if (t,i,j) in matchVars) <= 1 for t in range(T+K+1)), "only match with one other pair")
    model.addConstrs((quicksum(matchVars[t,i,j] for t in range(T+K+1) for j in range(K+1) if (t,i,j) in matchVars) + quicksum(matchVars[t,j,i] for t in range(T+K+1) for j in range(K+1) if (t,j,i) in matchVars) + quicksum(matchVars[i+T,j,jp] for j in range(K+1) \
                               for jp in range(K+1) if (i+T,j,jp) in matchVars) <= 1 for i in range(1,K+1)), "symmetry")
    
    if args.quality:
        obj = quicksum(matchVars[v]*matches[v] for v in matchVars)
    else:
        obj = quicksum(matchVars[v]*COUNT(v) for v in matchVars)
    model.setObjective(obj, GRB.MAXIMIZE) 
    model.optimize()
    for v in matchVars:
        if round(matchVars[v].X) != 0:
            count += COUNT(v)
            quality += matches[v]
            if v[2] == 0:
                agentInfo += "I" + str(v[0]-T) + "\t" + str(T+1) + "\t" + str(directed_matches[v[1]+T,v[0]]) + "\t" \
               + "I" + "\t" + str(directed_matches[v[0],v[1]+T]) + "\t" + "I" + "\t" + str(beta[v[0]-T]) + "\n"
                agentInfo += "I" + str(v[1]) + "\t" + str(T+1) + "\t" + str(directed_matches[v[0],v[1]+T]) + "\t" \
               + "I" + "\t" + str(directed_matches[v[1]+T,v[0]]) + "\t" + "I" + "\t" + str(beta[v[1]]) + "\n"
                del beta[v[0]-T]
                del beta[v[1]]
            else:
                agentInfo += "I" + str(v[0]-T) + "\t" + str(T+1) + "\t" + str(directed_matches[v[2]+T,v[0]]) + "\t" \
               + "I" + "\t" + str(directed_matches[v[0],v[1]+T]) + "\t" + "I" + "\t" + str(beta[v[0]-T]) + "\n"
                agentInfo += "I" + str(v[1]) + "\t" + str(T+1) + "\t" + str(directed_matches[v[0],v[1]+T]) + "\t" \
               + "I" + "\t" + str(directed_matches[v[1]+T,v[2]+T]) + "\t" + "I" + "\t" + str(beta[v[1]]) + "\n"
                agentInfo += "I" + str(v[2]) + "\t" + str(T+1) + "\t" + str(directed_matches[v[1]+T,v[2]+T]) + "\t" \
               + "I" + "\t" + str(directed_matches[v[2]+T,v[0]]) + "\t" + "I" + "\t" + str(beta[v[2]]) + "\n"
                del beta[v[0]-T]
                del beta[v[1]]
                del beta[v[2]]

               
            bt1 = getBloodTypes(demo[v[0]-1])
            bt2 = getBloodTypes(demo[v[1] + T - 1])
            graph += "edge [color="+graph_colors[bt1[1]] + "];\n"
            graph += "I" + str(v[0]-T) + " [color="+graph_colors[bt1[0]]+"];\n"
            graph += "I" + str(v[1]-1) + " [color="+graph_colors[bt2[0]]+"];\n"
            graph += "I" + str(v[0]-T) + " -> I" + str(v[1]-1) + ";\n"
            graph += "edge [color="+graph_colors[bt2[1]] + "];\n"
            graph += "node [color="+graph_colors[bt2[0]]+"];\n"
            graph += "I" + str(v[1]-1) + " -> I" + str(v[0]-T) + ";\n"
    del beta[0]
    for i in beta:
        agentInfo += "I" + str(i) + "\t" + str(T+2) + "\t" + str(0) + "\t" \
        + "N" + "\t" + str(0) + "\t" + "N" + "\t" + str(beta[i]) + "\n"

            
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
       
