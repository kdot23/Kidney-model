#!/usr/bin/env python2
# -*- coding: utf-8 -*-
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
dataIndex = 0

for fn in args.inputFiles:    
    with open(fn, 'rb') as f:
        d = pickle.load(f)
    num_incompat = d[0]
    num_compat = d[1]
    num_pairs = num_incompat + num_compat
    matches = d[2]
    T = num_compat
    K = num_incompat
    demo = d[4]
    
    quality = 0
    count = 0
    num_compat_to_self = 0
    num_compat_to_incompat = 0
    num_incompat_to_compat = 0
    num_incompat_to_incompat = 0
    
    graph = "digraph G { \n"
    #greedy algorithm for compatible pairs by order of index
    used_incompat = set()
    for i in range(1,num_compat+1):
        values = {j:matches[i,j] for j in range(K+1) if (i,j) in matches and j not in used_incompat}
        max_index = max(values, key=values.get)
        """
        for j in range(1,num_incompat+1):
            if matches[i,j] > matches[i][max_index] and j not in used_incompat and \
            sum(k>0 and k not in used_incompat for k in matches[i+T-1]) < args.n:
                max_index = j
        quality += matches[i][max_index]
        """
        quality += matches[i,max_index]
        count += COUNT((i,max_index))
        if max_index != 0: 
            used_incompat.add(max_index)
        if max_index == 0:
            bt = getBloodTypes(demo[i-1])
            graph += "edge [color="+graph_colors[bt[1]] + "];\n"
            graph += "node [color="+graph_colors[bt[0]]+"];\n"
            graph += "C" + str(i) + " -> C" + str(i) + ";\n"
            
        else:
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
            count += count(v)
            quality += matches[v]
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

if args.output:
    with open(args.output, 'w') as f:
        f.write(results)
else:
    print results
    
