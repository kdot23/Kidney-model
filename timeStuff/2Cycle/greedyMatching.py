#!/usr/bin/env python2
# -*- coding: utf-8 -*-
import pickle
import argparse
import numpy as np
from sets import Set
import os
import pulp

parser = argparse.ArgumentParser(description="Optimizes Kidney Exchange given by input file using a simple greedy algorithm")
parser.add_argument('--inputFiles', nargs='+', default = ["data.dat"], help="List of .dat files to be used as input. List of number of \
                    incompatible pairs, number of compatible pairs, list of quality(egs) of all possible pairs \
                    and demographic information. File created in KidneyDataGen")
parser.add_argument('--quality', action='store_true', help="Optimize for quality")
parser.add_argument('-o', '--output', help='write results to this file (.csv)')
parser.add_argument('--agents', help='output the quality of each agent to this file (.csv)')
parser.add_argument("-n", type = int, default = 2, help = "max number of connections incompatibles can be matched to and still removed from pool")
parser.add_argument("--graph", help = "output a graphviz representation")
parser.add_argument('-q', '--cadence', type=int)
parser.add_argument('--incompatible_online', type=float, help='threshhold for probability an incompatible stays before it is matched')
parser.add_argument('--gamma', default=.9, type=float, help='gamma value used for calculating survival')
args=parser.parse_args()

graph_colors = ["red", "blue", "green", "black"]

def COUNT(v):
    if v[1] == 0:
        return 2
    return 1

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
agentInfo = ''
dataIndex = 0

for fn in args.inputFiles:    
    with open(fn, 'rb') as f:
        d = pickle.load(f)
    I = d[0]
    C = d[1]
    num_pairs = I+C
    T = d[2]
    matches = d[3]
    demo = d[5]
    directed_matches = d[7]
    departure_times = d[8]

    available_incompat = set()
    unmatched_incompat = set()
    arriving_incompat = {}
    arriving_compat = {}
    
    quality = 0
    count = 0

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
        unmatched_incompat = unmatched_incompat.union(departing_incompat)
        if t in arriving_incompat:
            available_incompat = available_incompat.union(arriving_incompat[t])


        if t in arriving_compat:
            #do compatible matching stuff
            for i in arriving_compat[t]:
                values = {j:matches[i,j] for j in available_incompat.union(set([0])) if (i,j) in matches}
                max_index = max(values, key=values.get)
                if max_index == 0:
                    count += 1
                    agentInfo += "C" + str(i) + "\t" + str(t) + "\t" + str(directed_matches[i,0]) + "\t" \
                    + "C" + "\t" + str(directed_matches[i,0]) + "\t" + "C" + "\t" +  "\n"
                else:
                    available_incompat.remove(max_index)
                    count += 2
                    agentInfo += "C" + str(i) + "\t" + str(t) + "\t" + str(directed_matches[max_index+C,i]) + "\t" \
                    + "I" + "\t" + str(directed_matches[i,max_index+C]) + "\t" + "I" + "\t" +  "\n"
                    agentInfo += "I" + str(max_index) + "\t" + str(t) + "\t" + str(directed_matches[i,max_index+C]) + "\t" \
                    + "C" + "\t" + str(directed_matches[max_index+C,i]) + "\t" + "C" + "\t" +  "\n"
                quality += matches[i,max_index]


        if args.incompatible_online:
            probs = {i:demo[i+C-1][19]*args.gamma**num_rounds_present[i] for i in available_incompat}
            for i in probs:
                if i not in available_incompat: continue
                if probs[i] < args.incompatible_online:
                    values = {j:matches[i+C,j]  for j in available_incompat if (i+C,j) in matches}
                    if len(values) == 0: continue
                    max_index = max(values, key=values.get)
                    if values[max_index] > 0:
                        available_incompat.remove(i)
                        available_incompat.remove(max_index)
                        quality += matches[i+C,max_index]
                        count += COUNT((i+C,max_index))
                        agentInfo += "I" + str(i) + "\t" + str(t) + "\t" + str(directed_matches[max_index+C,i+C]) + "\t" \
                        + "I" + "\t" + str(directed_matches[i+C,max_index+C]) + "\t" + "I" + "\t" +  "\n"
                        agentInfo += "I" + str(max_index) + "\t" + str(t) + "\t" + str(directed_matches[i+C,max_index+C]) + "\t" \
                        + "I" + "\t" + str(directed_matches[max_index+C,i+C]) + "\t" + "I" + "\t" +  "\n"
        if args.cadence and  t%args.cadence==0:
            #do incompatible matching stuff
            matchVars = [(i+C,j) for i in available_incompat for j in available_incompat if (i+C, j) in matches]

            x = pulp.LpVariable.dicts('match', matchVars, lowBound = 0, upBound = 1, cat = pulp.LpInteger)
            matchVars = set(matchVars)
            model = pulp.LpProblem('incompatible matching', pulp.LpMaximize)
            if args.quality:
                model += lpSum(x[v]*matches[v] for v in matchVars)
            else:
                model += lpSum(x[v]*COUNT(v) for v in matchVars)
            for t in range(C+1,C+I+1):
                model += lpSum(x[t,i] for i in range(1,I+1) if (t,i) in matchVars) <= 1, 'compatible match with 1 for '+str(t) 
            for i in range(1,I+1):
                model += lpSum(x[t,i] for t in range(C+1,C+I+1) if (t,i) in matchVars) + lpSum(x[i+C,j] for j in range(1,I+1) if (i+C,j) in matchVars) <= 1, 'symetry '+str(i) 
            model.solve()
            count += sum(COUNT(v)*x[v] for v in matchVars)
            quality += sum(matches[v]*x[v] for v in matchVars)

            for v in matchVars:
                if round(x[v]) != 0:
                    available_incompat.remove(v[0]-C)
                    available_incompat.remove(v[1])
                    #Agent Info Stuff
                    agentInfo += "I" + str(v[0]-C) + "\t" + str(t) + "\t" + str(directed_matches[v[1]+C,v[0]]) + "\t" \
                    + "I" + "\t" + str(directed_matches[v[0],v[1]+C]) + "\t" + "I" + "\t" +  "\n"
                    agentInfo += "I" + str(v[1]) + "\t" + str(t) + "\t" + str(directed_matches[v[0],v[1]+C]) + "\t" \
                    + "I" + "\t" + str(directed_matches[v[1]+C,v[0]]) + "\t" + "I" + "\t" +  "\n"

    unmatched_incompat = unmatched_incompat.union(available_incompat)
    for i in unmatched_incompat:
        agentInfo += "I" + str(i) + "\t" + str(T) + "\t" + str(0) + "\t" \
        + "N" + "\t" + str(0) + "\t" + "N" + "\n"

    results += str(count) + '\t' + str(quality) + '\n'


"""
    graph = "digraph G { \n"
    #greedy algorithm for compatible pairs by order of index
    used_incompat = set()
    for i in range(1,num_compat+1):
        values = {j:matches[i,j] for j in range(K+1) if (i,j) in matches and j not in used_incompat}
        max_index = max(values, key=values.get)

        quality += matches[i,max_index]
        count += COUNT((i,max_index))
        if max_index != 0: 
            used_incompat.add(max_index)
        if max_index == 0:
            agentInfo += "C" + str(i) + "\t" + str(i) + "\t" + str(directed_matches[i,0]) + "\t" \
           + "C" + "\t" + str(directed_matches[i,0]) + "\t" + "C" + "\n"
            bt = getBloodTypes(demo[i-1])
            graph += "edge [color="+graph_colors[bt[1]] + "];\n"
            graph += "node [color="+graph_colors[bt[0]]+"];\n"
            graph += "C" + str(i) + " -> C" + str(i) + ";\n"
            
        else:
            agentInfo += "C" + str(i) + "\t" + str(i) + "\t" + str(directed_matches[max_index+T,i]) + "\t" \
           + "I" + "\t" + str(directed_matches[i,max_index+T]) + "\t" + "I" + "\n"
            agentInfo += "I" + str(max_index) + "\t" + str(i) + "\t" + str(directed_matches[i,max_index+T]) + "\t" \
           + "C" + "\t" + str(directed_matches[max_index+T,i]) + "\t" + "C" + "\n"
            bt1 = getBloodTypes(demo[i-1])
            bt2 = getBloodTypes(demo[max_index + T - 1])
            graph += "edge [color="+graph_colors[bt1[1]] + "];\n"
            graph += "C" + str(i) + " [color="+graph_colors[bt1[0]]+"];\n"
            graph += "I" + str(max_index) + " [color="+graph_colors[bt2[0]]+"];\n"
            graph += "C" + str(i) + " -> I" + str(max_index) + ";\n"
            graph += "edge [color="+graph_colors[bt2[1]] + "];\n"
            graph += "I" + str(max_index) + " -> C" + str(i) + ";\n"
    

    #gurobi optimization for remaining incompatible pairs
    model = Model('Kideny Optimizer')
    matchVars = {}
    for i in range(K):
        if i+1 in used_incompat: continue
        for j in range(K):
            if j+1 not in used_incompat and (i+T+1,j+1) in matches:
                matchVars[i+T+1,j+1] = model.addVar(vtype=GRB.BINARY, lb=0, ub=1, name=str((i+T+1,j+1)))
            
        
    
    model.addConstrs((quicksum(matchVars[i,j] for j in range(num_incompat+1) if (i,j) in matchVars) <= 1 \
                      for i in range(T+1,K+T+1)), "Only match with one pair")
    
    model.addConstrs((quicksum(matchVars[t,i] for t in range(1,T+K+1) if (t,i) in matchVars) + quicksum(matchVars[i+T,j] for j in range(1,K+1) \
            if (i+T,j) in matchVars) <= 1 for i in range(1,K+1)), "Symetry")

    if (args.quality):
        obj = quicksum(matchVars[i,j]*matches[i,j] for i in range(T+1,T+K+1) for j in range(1,K+1) if (i,j) in matchVars)
    else:
        obj = quicksum(COUNT((i,j))*matchVars[i,j] for i in range(T+1,T+K+1) for j in range(1,num_incompat+1) if (i,j) in matchVars)
        
    model.setObjective(obj, GRB.MAXIMIZE) 
    model.optimize()
    
    for v in matchVars:
        if round(matchVars[v].X) != 0:
            count += COUNT(v)
            quality += matches[v]
            agentInfo += "I" + str(v[0]-T) + "\t" + str(T+1) + "\t" + str(directed_matches[v[1]+T,v[0]]) + "\t" \
            + "N" + "\t" + str(directed_matches[v[0],v[1]+T]) + "\t" + "I" + "\n"
            agentInfo += "I" + str(v[1]) + "\t" + str(T+1) + "\t" + str(directed_matches[v[0],v[1]+T]) + "\t" \
            + "N" + "\t" + str(directed_matches[v[1]+T,v[0]]) + "\t" + "I" + "\n"
            used_incompat.add(v[0]-T)
            used_incompat.add(v[1])
            bt1 = getBloodTypes(demo[v[0]-1])
            bt2 = getBloodTypes(demo[v[1]+T-1])
            graph += "edge [color="+graph_colors[bt1[1]] + "];\n"
            graph += "I" + str(v[0]) + " [color="+graph_colors[bt1[0]]+"];\n"
            graph += "I" + str(v[1]) + " [color="+graph_colors[bt2[0]]+"];\n"
            graph += "I" + str(v[0]) + " -> I" + str(v[1]) + ";\n"
            graph += "edge [color="+graph_colors[bt2[1]] + "];\n"
            graph += "I" + str(v[1]) + " -> I" + str(v[0]) + ";\n"
    graph += "}"
    if args.graph:
        with open(args.graph+str(dataIndex)+".gv", 'w') as f:
            f.write(graph)
        os.system('dot -Tpdf ' + args.graph + str(dataIndex) + ".gv -o " + args.graph + str(dataIndex) +".pdf")
            
    dataIndex += 1
    results += str(count) + "\t" + str(quality) + "\n"
    for i in range(1,num_incompat+1):
        if i not in used_incompat:
            agentInfo += "I" + str(i) + "\t" + str(T+2) + "\t" + str(0) + "\t" \
            + "N" + "\t" + str(0) + "\t" + "N" + "\n"
            """

if args.output:
    with open(args.output, 'w') as f:
        f.write(results)
else:
    print results
    
if args.agents:
    with open(args.agents, 'w') as f:
        f.write(agentInfo)
