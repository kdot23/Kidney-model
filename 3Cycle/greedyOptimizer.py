#!/usr/bin/env python2
# -*- coding: utf-8 -*-
import json
import pickle
import argparse
from gurobipy import *
import numpy as np
from sets import Set
import os

parser = argparse.ArgumentParser(description="Optimizes Kidney Exchange given by input file using a simple greedy algorithm")
parser.add_argument('--inputFiles', nargs='+', default = ["data.dat"], help="List of .dat files to be used as input. List of number of \
                    incompatible pairs, number of compatible pairs, list of quality(egs) of all possible pairs \
                    and demographic information. File created in KidneyDataGen")
parser.add_argument('--quality', action='store_true', help="Optimize for quality")
parser.add_argument('-o', '--output', help='write results to this file (.csv)')
parser.add_argument("-n", type = int, default = 2, help = "max number of connections incompatibles can be matched to and still removed from pool")
parser.add_argument("--graph", help = "output a graphviz representation")
args=parser.parse_args()

graph_colors = ["red", "blue", "green", "black"]

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


data = []
for fn in args.inputFiles:
    with open (fn, 'rb') as f:
        data.append(pickle.load(f))

results = ''
dataIndex = 0

for d in data:    
    num_incompat = d[0]
    num_compat = d[1]
    num_pairs = num_incompat + num_compat
    matches = d[2]
    T = num_compat
    demo = d[3]
    
    quality = 0
    num_compat_to_self = 0
    num_compat_to_incompat = 0
    num_incompat_to_compat = 0
    num_incompat_to_incompat = 0
    
    graph = "digraph G { \n"
    #greedy algorithm for compatible pairs by order of index
    used_incompat = set()
    for i in range(num_compat):
        max_index = max((v for v in matches if v[0] == i and v[1] not in used_incompat and v[2] not in used_incompat), key=matches.get)
        quality += matches[max_index]
        if max_index[1] != 0:
            used_incompat.add(max_index[1])
        if max_index[2] != 0:
            used_incompat.add(max_index[2])
        if max_index[1] == 0:
            num_compat_to_self += 1
            bt = getBloodTypes(demo[i])
            graph += "edge [color="+graph_colors[bt[1]] + "];\n"
            graph += "node [color="+graph_colors[bt[0]]+"];\n"
            graph += "C" + str(i) + " -> C" + str(i) + ";\n"
            
        elif max_index[2] == 0:
            num_compat_to_incompat += 1
            num_incompat_to_compat += 1
            bt1 = getBloodTypes(demo[i])
            bt2 = getBloodTypes(demo[max_index[1] + T - 1])
            graph += "edge [color="+graph_colors[bt1[1]] + "];\n"
            graph += "C" + str(i) + " [color="+graph_colors[bt1[0]]+"];\n"
            graph += "I" + str(max_index[1]-1) + " [color="+graph_colors[bt2[0]]+"];\n"
            graph += "C" + str(i) + " -> I" + str(max_index[1]-1) + ";\n"
            graph += "edge [color="+graph_colors[bt2[1]] + "];\n"
            graph += "I" + str(max_index[1]-1) + " -> C" + str(i) + ";\n"
        else:
            num_compat_to_incompat += 1
            num_incompat_to_incompat += 1
            num_incompat_to_compat += 1
    
    #gurobi optimization for remaining incompatible pairs
    model = Model('Kideny Optimizer')
    matchVars = {v:model.addVar(vtype = GRB.CONTINUOUS, lb = 0, ub=1,  name = "incompat_match_" + str(v)) for v in matches if v[0] >= T and \
            v[0]-T+1 not in used_incompat and v[1] not in used_incompat and v[2] not in used_incompat}

    model.addConstrs((quicksum(matchVars[t,i,j] for i in range(num_incompat+1) for j in range(num_incompat+1) if (t,i,j) in matchVars) <= 1 \
                      for t in range(num_pairs)), "Only match with one pair")
    
    model.addConstrs((quicksum(matchVars[t,i,j] for t in range(num_pairs) for j in range(num_incompat+1) if (t,i,j) in matchVars) \
                     + quicksum(matchVars[t,j,i] for t in range(num_pairs) for j in range(1, num_incompat+1) if (t,j,i) in matchVars)\
                     + quicksum(matchVars[i+T-1,j,jp] for j in range(1,num_incompat+1) for jp in range(num_incompat+1) if (i+T-1,j,jp) in matchVars)\
                     <= 1 for i in range(1,num_incompat+1)), "undirected graph")

    if (args.quality):
        obj = quicksum(matchVars[v]*matches[v] for v in matchVars)
    else:
        obj = quicksum(COUNT(v)*matchVars[v] for v in matchVars)


    model.setObjective(obj, GRB.MAXIMIZE) 
    model.optimize()
    
    quality += obj.getValue()/2
    for v in matchVars:
        if matchVars[v].X != 0:
            num_incompat_to_incompat += 1
            bt1 = getBloodTypes(demo[v[0] + T])
            bt2 = getBloodTypes(demo[v[1] + T])
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
            
    num_matches = num_compat_to_self + num_compat_to_incompat + num_incompat_to_compat + num_incompat_to_incompat
    dataIndex += 1
    
    results += str(num_matches) + "\t" + str(quality) + "\n"

if args.output:
    with open(args.output, 'w') as f:
        f.write(results)
else:
    print results
    
