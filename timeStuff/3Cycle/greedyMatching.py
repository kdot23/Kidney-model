#!/usr/bin/env python2
# -*- coding: utf-8 -*-
import json
import pickle
import argparse
import numpy as np
from sets import Set
import os
import pulp
from pulp import lpSum

parser = argparse.ArgumentParser(description="Optimizes Kidney Exchange given by input file using a simple greedy algorithm")
parser.add_argument('--inputFiles', nargs='+', default = ["data.dat"], help="List of .dat files to be used as input. List of number of \
                    incompatible pairs, number of compatible pairs, list of quality(egs) of all possible pairs \
                    and demographic information. File created in KidneyDataGen")
parser.add_argument('--quality', action='store_true', help="Optimize for quality")
parser.add_argument('-o', '--output', help='write results to this file (.csv)')
parser.add_argument('--agents', help='output the quality of each agent to this file (.csv)')
parser.add_argument("-n", type = int, default = 2, help = "max number of connections incompatibles can be matched to and still removed from pool")
parser.add_argument("--graph", help = "output a graphviz representation")
parser.add_argument('-q', '--cadence',  type=int)
parser.add_argument('--incompatible_online', type=float, help='threshhold for probability an incompatible stays before it is matched')
parser.add_argument('--gamma', default=.9, type=float, help='gamma value used for calculating survival')
args=parser.parse_args()

graph_colors = ["red", "blue", "green", "black"]

def COUNT(v):
    if v[2] != 0:
        return 3
    if v[1] != 0:
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
agentInfo = '' #Pair id, time of match, get quality, get donor type, give quality, give recipient type
dataIndex = 0

for fn in args.inputFiles:
    with open(fn, 'rb') as f:
        d = pickle.load(f)
    I = d[0]
    C = d[1]
    num_pairs = I + C
    matches = d[4]
    directed_matches = d[7]
    T = d[2]
    demo = d[5]
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
            #Do compatible matching stuff
            for i in arriving_compat[t]:
                values = {(i,j,k):matches[i,j,k] for j in available_incompat.union(set([0])) for k in available_incompat.union(set([0])) 
                        if (i,j,k) in matches}
                max_index = max(values, key=values.get)
                if max_index[1] == 0:
                    count += 1
                    agentInfo += "C" + str(i) + "\t" + str(t) + "\t" + str(directed_matches[max_index[0],0]) + "\t" \
                    + "C" + "\t" + str(directed_matches[max_index[0],0]) + "\t" + "C" + "\n"
                elif max_index[2] == 0:
                    available_incompat.remove(max_index[1])
                    count += 2
                    agentInfo += "C" + str(i) + "\t" + str(t) + "\t" + str(directed_matches[max_index[1]+C,max_index[0]]) + "\t" \
                    + "I" + "\t" + str(directed_matches[max_index[0],max_index[1]+C]) + "\t" + "I" + "\n"
                    agentInfo += "I" + str(max_index[1]) + "\t" + str(t) + "\t" + str(directed_matches[max_index[0],max_index[1]+C]) + "\t" \
                    + "C" + "\t" + str(directed_matches[max_index[1]+C,max_index[0]]) + "\t" + "C" + "\n"
                else:
                    available_incompat.remove(max_index[1])
                    available_incompat.remove(max_index[2])
                    count += 3
                    agentInfo += "C" + str(i) + "\t" + str(t) + "\t" + str(directed_matches[max_index[2]+C,max_index[0]]) + "\t" \
                    + "I" + "\t" + str(directed_matches[max_index[0],max_index[1]+C]) + "\t" + "I" + "\n"
                    agentInfo += "I" + str(max_index[1]) + "\t" + str(t) + "\t" + str(directed_matches[max_index[0],max_index[1]+C]) + "\t" \
                    + "C" + "\t" + str(directed_matches[max_index[1]+C,max_index[2]+C]) + "\t" + "I" + "\n"
                    agentInfo += "I" + str(max_index[2]) + "\t" + str(t) + "\t" + str(directed_matches[max_index[1]+C,max_index[2]+C]) + "\t" \
                    + "I" + "\t" + str(directed_matches[max_index[2]+C,max_index[0]]) + "\t" + "C" + "\n"
                quality += matches[max_index]

        if args.incompatible_online:
            probs = {i:demo[i+C-1][19]*args.gamma**num_rounds_present[i] for i in available_incompat}
            for i in probs:
                if i not in available_incompat: continue
                if probs[i] < args.incompatible_online:
                    values = {(i+C,j,k):matches[i+C,j, k] for j in available_incompat for k in available_incompat if (i+C,j,k) in matches}
                    if len(values) == 0: continue
                    max_index = max(values, key=values.get)
                    if values[max_index] > 0:
                        available_incompat.remove(i)
                        available_incompat.remove(max_index[1])
                        quality += matches[max_index]
                        count += COUNT(max_index)
                        if max_index[2] == 0:
                            agentInfo += "I" + str(i) + "\t" + str(t) + "\t" + str(directed_matches[max_index[1]+C,max_index[0]]) + "\t" \
                            + "I" + "\t" + str(directed_matches[max_index[0],max_index[1]+C]) + "\t" + "I" + "\n"
                            agentInfo += "I" + str(max_index[1]) + "\t" + str(t) + "\t" + str(directed_matches[max_index[0],max_index[1]+C]) + "\t" \
                            + "I" + "\t" + str(directed_matches[max_index[1]+C,max_index[0]]) + "\t" + "I" + "\n"
                        else:
                            available_incompat.remove(max_index[2])
                            agentInfo += "I" + str(i) + "\t" + str(t) + "\t" + str(directed_matches[max_index[2]+C,max_index[0]]) + "\t" \
                            + "I" + "\t" + str(directed_matches[max_index[0],max_index[1]+C]) + "\t" + "I" + "\n"
                            agentInfo += "I" + str(max_index[1]) + "\t" + str(t) + "\t" + str(directed_matches[max_index[0],max_index[1]+C]) + "\t" \
                            + "I" + "\t" + str(directed_matches[max_index[1]+C,max_index[2]+C]) + "\t" + "I" + "\n"
                            agentInfo += "I" + str(max_index[2]) + "\t" + str(t) + "\t" + str(directed_matches[max_index[1]+C,max_index[2]+C]) + "\t" \
                            + "I" + "\t" + str(directed_matches[max_index[2]+C,max_index[0]]) + "\t" + "I" + "\n"
        if args.cadence and t%args.cadence==0:
            #Do incompatible matching stuff
            model = pulp.LpProblem('incompatible matching', pulp.LpMaximize)
            matchVars = [(i+C,j,k) for i in available_incompat for j in available_incompat for k in available_incompat \
                    if (i+C,j,k) in matches]
            x = pulp.LpVariable.dicts('match',matchVars,lowBound = 0, upBound = 1, cat = pulp.LpInteger)
            matchVars = set(matchVars)
            if args.quality:
                model += lpSum(x[v]*matches[v] for v in matchVars)
            else:
                model += lpSum(x[v]*COUNT(v) for v in matchVars)
            for t in available_incompat:
                model += lpSum(x[t+C,i,j] for i in available_incompat for j in available_incompat if (t+C,i,j) in matchVars) <= 1,\
                        'match with one '+str(t)
            for i in available_incompat:
                model += lpSum(x[t+C,i,j] for t in available_incompat for j in available_incompat if (t+C,i,j) in matchVars) + \
                        lpSum(x[t+C,j,i] for t in available_incompat for j in available_incompat if (t+C,j,i) in matchVars) + \
                        lpSum(x[i+C,j,k] for j in available_incompat for k in available_incompat if (i+C,j,k) in matchVars) <= 1, \
                        'symetry ' + str(i)

            model.solve()
            count += sum(COUNT(v)*x[v].value() for v in matchVars)
            quality += sum(matches[v]*x[v].value() for v in matchVars)

            for v in matchVars:
                if round(x[v].value()) != 0:
                    available_incompat.remove(v[0]-C)
                    available_incompat.remove(v[1])
                    if v[2] == 0:
                        agentInfo += "I" + str(v[0]-C) + "\t" + str(t) + "\t" + str(directed_matches[v[1]+C,v[0]]) + "\t" \
                        + "I" + "\t" + str(directed_matches[v[0],v[1]+C]) + "\t" + "I" + "\n"
                        agentInfo += "I" + str(v[1]) + "\t" + str(t) + "\t" + str(directed_matches[v[0],v[1]+C]) + "\t" \
                        + "I" + "\t" + str(directed_matches[v[1]+C,v[0]]) + "\t" + "I" + "\n"
                    else:
                        available_incompat.remove(v[2])
                        agentInfo += "I" + str(v[0]-C) + "\t" + str(t) + "\t" + str(directed_matches[v[2]+C,v[0]]) + "\t" \
                        + "I" + "\t" + str(directed_matches[v[0],v[1]+C]) + "\t" + "I" + "\n"
                        agentInfo += "I" + str(v[1]) + "\t" + str(t) + "\t" + str(directed_matches[v[0],v[1]+C]) + "\t" \
                        + "I" + "\t" + str(directed_matches[v[1]+C,v[2]+C]) + "\t" + "I" + "\n"
                        agentInfo += "I" + str(v[2]) + "\t" + str(t) + "\t" + str(directed_matches[v[1]+C,v[2]+C]) + "\t" \
                        + "I" + "\t" + str(directed_matches[v[2]+C,v[0]]) + "\t" + "I" + "\n"

    unmatched_incompat = unmatched_incompat.union(available_incompat)
    results += str(count) + '\t' + str(quality) + '\n'
    for i in unmatched_incompat:
        agentInfo += "I" + str(i) + "\t" + str(T) + "\t" + str(0) + "\t" \
        + "N" + "\t" + str(0) + "\t" + "N" + "\n"


    """
    
    graph = "digraph G { \n"
    #greedy algorithm for compatible pairs by order of index
    used_incompat = set()
    for i in range(1,num_compat+1):
        max_index = max((v for v in matches if v[0] == i and v[1] not in used_incompat and v[2] not in used_incompat), key=matches.get)
        quality += matches[max_index]

        count += COUNT(max_index)
        if max_index[1] != 0:
            used_incompat.add(max_index[1])
        if max_index[2] != 0:
            used_incompat.add(max_index[2])
        
        if max_index[1] == 0:
            agentInfo += "C" + str(i) + "\t" + str(i) + "\t" + str(directed_matches[max_index[0],0]) + "\t" \
           + "C" + "\t" + str(directed_matches[max_index[0],0]) + "\t" + "C" + "\n"
            bt = getBloodTypes(demo[i])
            graph += "edge [color="+graph_colors[bt[1]] + "];\n"
            graph += "node [color="+graph_colors[bt[0]]+"];\n"
            graph += "C" + str(i) + " -> C" + str(i) + ";\n"
            
        elif max_index[2] == 0:
            agentInfo += "C" + str(i) + "\t" + str(i) + "\t" + str(directed_matches[max_index[1]+T,max_index[0]]) + "\t" \
           + "I" + "\t" + str(directed_matches[max_index[0],max_index[1]+T]) + "\t" + "I" + "\n"
            agentInfo += "I" + str(max_index[1]) + "\t" + str(i) + "\t" + str(directed_matches[max_index[0],max_index[1]+T]) + "\t" \
           + "C" + "\t" + str(directed_matches[max_index[1]+T,max_index[0]]) + "\t" + "C" + "\n"
           
            bt1 = getBloodTypes(demo[i])
            bt2 = getBloodTypes(demo[max_index[1] + T - 1])
            graph += "edge [color="+graph_colors[bt1[1]] + "];\n"
            graph += "C" + str(i) + " [color="+graph_colors[bt1[0]]+"];\n"
            graph += "I" + str(max_index[1]-1) + " [color="+graph_colors[bt2[0]]+"];\n"
            graph += "C" + str(i) + " -> I" + str(max_index[1]-1) + ";\n"
            graph += "edge [color="+graph_colors[bt2[1]] + "];\n"
            graph += "I" + str(max_index[1]-1) + " -> C" + str(i) + ";\n"

        else:
            agentInfo += "C" + str(i) + "\t" + str(i) + "\t" + str(directed_matches[max_index[2]+T,max_index[0]]) + "\t" \
           + "I" + "\t" + str(directed_matches[max_index[0],max_index[1]+T]) + "\t" + "I" + "\n"
            agentInfo += "I" + str(max_index[1]) + "\t" + str(i) + "\t" + str(directed_matches[max_index[0],max_index[1]+T]) + "\t" \
           + "C" + "\t" + str(directed_matches[max_index[1]+T,max_index[2]+T]) + "\t" + "I" + "\n"
            agentInfo += "I" + str(max_index[2]) + "\t" + str(i) + "\t" + str(directed_matches[max_index[1]+T,max_index[2]+T]) + "\t" \
           + "I" + "\t" + str(directed_matches[max_index[2]+T,max_index[0]]) + "\t" + "C" + "\n"
        
    #gurobi optimization for remaining incompatible pairs
    model = Model('Kideny Optimizer')
    matchVars = {v:model.addVar(vtype = GRB.BINARY, lb = 0, ub=1,  name = "incompat_match_" + str(v)) for v in matches if v[0] > T and \
            v[0]-T not in used_incompat and v[1] not in used_incompat and v[2] not in used_incompat}

    model.addConstrs((quicksum(matchVars[t,i,j] for i in range(num_incompat+1) for j in range(num_incompat+1) if (t,i,j) in matchVars) <= 1 \
                      for t in range(1,num_pairs+1)), "Only match with one pair")
    
    model.addConstrs((quicksum(matchVars[t,i,j] for t in range(1,num_pairs+1) for j in range(num_incompat+1) if (t,i,j) in matchVars) \
                     + quicksum(matchVars[t,j,i] for t in range(1,num_pairs+1) for j in range(1, num_incompat+1) if (t,j,i) in matchVars)\
                     + quicksum(matchVars[i+T,j,jp] for j in range(1,num_incompat+1) for jp in range(num_incompat+1) if (i+T,j,jp) in matchVars)\
                     <= 1 for i in range(1,num_incompat+1)), "undirected graph")

    if (args.quality):
        obj = quicksum(matchVars[v]*matches[v] for v in matchVars)
    else:
        obj = quicksum(COUNT(v)*matchVars[v] for v in matchVars)


    model.setObjective(obj, GRB.MAXIMIZE) 
    model.optimize()
    

    for v in matchVars:
        if round(matchVars[v].X) != 0:
            quality += matches[v]
            count += COUNT(v)
            used_incompat.add(v[0]-T)
            used_incompat.add(v[1])
            if v[2] == 0:
                 agentInfo += "I" + str(v[0]-T) + "\t" + str(T+1) + "\t" + str(directed_matches[v[1]+T,v[0]]) + "\t" \
               + "I" + "\t" + str(directed_matches[v[0],v[1]+T]) + "\t" + "I" + "\n"
                 agentInfo += "I" + str(v[1]) + "\t" + str(T+1) + "\t" + str(directed_matches[v[0],v[1]+T]) + "\t" \
               + "I" + "\t" + str(directed_matches[v[1]+T,v[0]]) + "\t" + "I" + "\n"
            else:
                 used_incompat.add(v[2])
                 agentInfo += "I" + str(v[0]-T) + "\t" + str(T+1) + "\t" + str(directed_matches[v[2]+T,v[0]]) + "\t" \
               + "I" + "\t" + str(directed_matches[v[0],v[1]+T]) + "\t" + "I" + "\n"
                 agentInfo += "I" + str(v[1]) + "\t" + str(T+1) + "\t" + str(directed_matches[v[0],v[1]+T]) + "\t" \
               + "I" + "\t" + str(directed_matches[v[1]+T,v[2]+T]) + "\t" + "I" + "\n"
                 agentInfo += "I" + str(v[2]) + "\t" + str(T+1) + "\t" + str(directed_matches[v[1]+T,v[2]+T]) + "\t" \
               + "I" + "\t" + str(directed_matches[v[2]+T,v[0]]) + "\t" + "I" + "\n"
               
            bt1 = getBloodTypes(demo[v[0]-1])
            bt2 = getBloodTypes(demo[v[1] + T - 1])
            graph += "edge [color="+graph_colors[bt1[1]] + "];\n"
            graph += "I" + str(v[0]) + " [color="+graph_colors[bt1[0]]+"];\n"
            graph += "I" + str(v[1]) + " [color="+graph_colors[bt2[0]]+"];\n"
            graph += "I" + str(v[0]) + " -> I" + str(v[1]) + ";\n"
            graph += "edge [color="+graph_colors[bt2[1]] + "];\n"
            graph += "I" + str(v[1]) + " -> I" + str(v[0]) + ";\n"

    for i in range(1,num_incompat+1):
        if i not in used_incompat:
             agentInfo += "I" + str(i) + "\t" + str(T+2) + "\t" + str(0) + "\t" \
        + "N" + "\t" + str(0) + "\t" + "N" + "\n"

    graph += "}"
    if args.graph:
        with open(args.graph+str(dataIndex)+".gv", 'w') as f:
            f.write(graph)
        os.system('dot -Tpdf ' + args.graph + str(dataIndex) + ".gv -o " + args.graph + str(dataIndex) +".pdf")
            
    dataIndex += 1
    
    results += str(count) + "\t" + str(quality) + "\n"
    """

if args.output:
    with open(args.output, 'w') as f:
        f.write(results)
        
else:
    print results

if args.agents:
    with open(args.agents, 'w') as f:
        f.write(agentInfo)
    
