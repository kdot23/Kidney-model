#!/usr/bin/env python2
# -*- coding: utf-8 -*-
import json
import argparse
from gurobipy import *
import numpy as np
from sets import Set
import os

parser = argparse.ArgumentParser(description="Optimizes Kidney Exchange given by input file using a simple greedy algorithm")
parser.add_argument('--inputFile', nargs='?', help="JSON File to be used as input. List of number of \
                    incompatible pairs, number of compatible pairs, list of quality(egs) of all possible pairs \
                    and demographic information. File created in KidneyDataGen")
parser.add_argument('--quality', action='store_true', help="Optimize for quality")
parser.add_argument('-i', '--inputDir', nargs='?', default='data', help='input directory to look for data files')
parser.add_argument('-o', '--output', help='write results to this file')
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

data = []
if args.inputFile:
    with open (args.inputFile, 'r') as f:
        data.append(json.load(f))
else:
    for fn in os.listdir(args.inputDir):
        if fn == '.DS_Store': continue
        with open(args.inputDir+'/'+fn, 'r') as f:
            data.append(json.load(f))

pastData = []

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
        max_index = 0
        for j in range(1,num_incompat+1):
            if matches[i][j] > matches[i][max_index] and j not in used_incompat and \
            sum(k>0 and k not in used_incompat for k in matches[i+T-1]) < args.n:
                max_index = j
        quality += matches[i][max_index]
        used_incompat.add(max_index)
        if max_index == 0:
            num_compat_to_self += 1
            bt = getBloodTypes(demo[i])
            graph += "edge [color="+graph_colors[bt[1]] + "];\n"
            graph += "node [color="+graph_colors[bt[0]]+"];\n"
            graph += "C" + str(i) + " -> C" + str(i) + ";\n"
            
        else:
            num_compat_to_incompat += 1
            num_incompat_to_compat += 1
            bt1 = getBloodTypes(demo[i])
            bt2 = getBloodTypes(demo[max_index + T - 1])
            graph += "edge [color="+graph_colors[bt1[1]] + "];\n"
            graph += "node [color="+graph_colors[bt1[0]]+"];\n"
            graph += "C" + str(i) + " -> I" + str(max_index-1) + ";\n"
            graph += "edge [color="+graph_colors[bt2[1]] + "];\n"
            graph += "node [color="+graph_colors[bt2[0]]+"];\n"
            graph += "I" + str(max_index-1) + " -> C" + str(i) + ";\n"
    
    #gurobi optimization for remaining incompatible pairs
    model = Model('Kideny Optimizer')
    matchVars = {}
    incompat_matches = []
    for i in range(num_incompat):
        incompat_matches.append(matches[i + num_compat][1:])

    for i in range(num_incompat):
        if i+1 in used_incompat: continue
        for j in range(num_incompat):
            if j+1 not in used_incompat and incompat_matches[i][j] != 0:
                matchVars[(i,j)] = model.addVar(vtype = GRB.CONTINUOUS, lb = 0, ub=1,  name = "incompat_match_" + str((i,j)))
    
    model.addConstrs((quicksum(matchVars[i,j] for j in range(num_incompat) if (i,j) in matchVars) <= 1 \
                      for i in range(num_incompat)), "Only match with one pair")
    
    model.addConstrs((matchVars[(i,j)] == matchVars[(j,i)] for i in range(num_incompat) for j in range(num_incompat) \
                                if (i,j) in matchVars and (j,i) in matchVars), "undirected graph")

    if (args.quality):
        obj = quicksum(matchVars[i,j]*matches[i][j] for i in range(num_incompat) for j in range(num_incompat) if (i,j) in matchVars and (j,i) in matchVars)
    else:
        obj = quicksum(matchVars[i,j] for i in range(num_incompat) for j in range(num_incompat+1) if (i,j) in matchVars and (j,i) in matchVars)
        
    model.setObjective(obj, GRB.MAXIMIZE) 
    model.optimize()
    
    quality += obj.getValue()
    for v in matchVars:
        if matchVars[v].X != 0:
            num_incompat_to_incompat += 1
            bt1 = getBloodTypes(demo[v[0] + T])
            bt2 = getBloodTypes(demo[v[1] + T])
            graph += "edge [color="+graph_colors[bt1[1]] + "];\n"
            graph += "node [color="+graph_colors[bt1[0]]+"];\n"
            graph += "I" + str(v[0]) + " -> I" + str(v[1]) + ";\n"
            graph += "edge [color="+graph_colors[bt2[1]] + "];\n"
            graph += "node [color="+graph_colors[bt2[0]]+"];\n"
            graph += "I" + str(v[1]) + " -> I" + str(v[0]) + ";\n"
            
    num_matches = num_compat_to_self + num_compat_to_incompat + num_incompat_to_compat + num_incompat_to_incompat
"""
    pastData.append((quality, num_matches))
    if args.output:
        with open(args.output, 'a') as f:
            f.write(str(quality) + "\t" + str(num_matches) + "\t" +  str(num_incompat_to_incompat) +"\n")
                     #num_compat_to_self, num_compat_to_incompat, num_incompat_to_compat, num_incompat_to_incompat))
print str(num_matches) + "\t" + str(quality) + "\n"
avgs = np.mean(pastData, axis=0)
stdevs = np.std(pastData, axis=0)
s = ''
results = ''
for std in stdevs:
    s += str(std)+"\t"
s += "\n"
results = s+"\n\n\n"+results
s = ''
for a in avgs:
    s += str(a)+"\t"
s+="\n"
results = s+results
"""

if args.output:
    with open(args.output, 'w') as f:
        f.write(str(num_matches) + "\t" + str(quality) + "\n")
else:
    print str(num_matches) + "\t" + str(quality) + "\n"
    
graph += "}"
if args.graph:
    with open(args.graph+".gv", 'w') as f:
        f.write(graph)
    os.system('dot -Tpdf ' + args.graph + ".gv -o " + args.graph +".pdf")
