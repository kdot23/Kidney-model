#inport modules
import numpy as np
import os
import json
from gurobipy import *
import argparse

parser = argparse.ArgumentParser(description="Optimizes Kidney Exchange given by input file")
parser.add_argument('--inputFile', nargs='?', help="JSON File to be used")
parser.add_argument('-i', '--inputDir', nargs='?', default='data', help='input directory to look for data files')
parser.add_argument('-o', '--outputFile', nargs='?', help="CSV File for results to be saved to")
parser.add_argument('--quality', action='store_true', help="Optimize for quality")
parser.add_argument('--count', action='store_true', help="Optimize for count")
parser.add_argument('-a', '--all', action='store_true', help="Perform all optimizations")

args=parser.parse_args()
args.count = args.count or args.all
args.quality = not args.count or args.quality or args.all
results = "#nm\tnum_C\tnum_CI\tnum_II\tq\tc_q\tci_q\tii_q\n"
if args.inputFile:
    with open(args.inputFile,'r') as f:
        data = [json.load(f)]
else:
    data = []
    for fn in os.listdir(args.inputDir):
        with open(args.inputDir+'/'+fn, 'r') as f:
            data.append(json.load(f))

pastData = []
for d in data:
    compat = d[0]
    incompat = d[1]
    Igraph = d[2]
    CIgraph = d[3]
    T = len(compat)
    
    model = Model('Kideny Optimizer')
    iMatches = {}
    ciMatches = {}
    #only create variables for edges which exist
    for i in range(T):
        for j in range(T):
            if Igraph[i][j] and Igraph[j][i]:
                iMatches[(i,j)] = model.addVar(vtype = GRB.CONTINUOUS, lb = 0, ub = 1, name = "I" + str((i,j)))
            if CIgraph[i][j] and incompat[j][2] > compat[i][2]:
                ciMatches[(i,j)] = model.addVar(vtype = GRB.CONTINUOUS, lb = 0, ub = 1, name = "CI" + str((i,j)))
        ciMatches[(i,-1)] = model.addVar(vtype = GRB.CONTINUOUS, lb = 0, ub = 1, name = "CI" + str((i,-1)))
    
    #first entry is which compatible pair, second entry is incompatible pair or -1 for self matching
    #ciMatches are compatible to incompatible or compatible to themselves (-1)
    #ciMatches = model.addVars(range(T), range(-1, T), vtype = GRB.BINARY, name="CIMatches")
    #iMatches are incompatible to incompatible
    #iMatches = model.addVars(range(T), range(T), vtype = GRB.BINARY, name = "IMatches")
    
    #compatible pair can only match with up to one incompatible pair or themselves
    model.addConstrs((quicksum(ciMatches[i,j] for j in range(-1,T) if (i,j) in ciMatches)  <= 1 for i in range(T)), "Compat  Matches")
    #incompatible pair can only match with up to one incompatible pair or one compatible pair
    model.addConstrs((quicksum(iMatches[i,j] for j in range(T) if (i,j) in iMatches) + quicksum(ciMatches[j,i] for j in range(T) if (j,i) in ciMatches) <= 1 for i in range(T)), "Incompat Matches")
    
    model.addConstrs((iMatches[i,j] == iMatches[j,i] for i in range(T) for j in range(T) if (i,j) in iMatches), "Symmetry")
    
    #compatible will only exchange if they have a higher quality score due to the exchange
    #model.addConstrs((ciMatches[i,j]*incompat[j][2] >= ciMatches[i,j]*compat[i][2] for i in range(T) for j in range(T)), " Compat Selfishness")
    
    #an edge between C and I can only exist if the pairs are compatible
    #model.addConstrs((CIgraph[i][j] >= ciMatches[i,j] for i in range(T) for j in range(T)), "Compatible to Incompatible Compatibility")
    
    #an edge between two Is can only exist if the pairs are compatible 
    #model.addConstrs((Igraph[i][j] >= iMatches[i,j] for i in range(T) for j in range(T)), "Incompatible Compatibility")
    
    #substitute this line for line 24 to allow for arbitrary length cycles in incompatible pool
    #model.addConstrs((quicksum(iMatches[i,j] for j in range(T) if (i,j) in iMatches) == quicksum(iMatches[j,i] for j in range(T)) for i in range(T)), "Give one to get one")
    
    #obj = quicksum(ciMatches[i,j] for i in range(T) for j in range(T) if (i,j) in ciMatches) + quicksum(ciMatches[i,-1] for i in range(T) ) + quicksum(iMatches[i,j] for i in range(T) for j in range(T) if (i,j) in iMatches)
    obj = quicksum((compat[i][2]+incompat[j][2])*ciMatches[i,j] for i in range(T) for j in range(T) if (i,j) in ciMatches) + quicksum(compat[i][2]*ciMatches[i,-1] for i in range(T) ) + quicksum(incompat[i][2]*iMatches[i,j] for i in range(T) for j in range(T) if (i,j) in iMatches)
    model.setObjective(obj, GRB.MAXIMIZE) 
    model.optimize()
    num_matches = 0
    quality = 0.
    c_quality = 0.
    ci_quality = 0.
    ii_quality = 0.
    num_II = 0
    num_CI = 0
    num_C = 0
    for v in model.getVars():
        if v.X != 0:
            #print str(v.Varname) + ' ' + str(v.X)
            num_matches += 1
    
    for var in iMatches:
        if (iMatches[var].X != 0):
            ii_quality += incompat[var[0]][2]
            num_II += 1
    for var in ciMatches:
        if (ciMatches[var].X != 0):
            if var[1] == -1:
                num_C += 1
                c_quality += compat[var[0]][2]
            else:
                num_CI += 1
                ci_quality += compat[var[0]][2] + incompat[var[1]][2]
    
    print num_matches
    print num_C
    print num_CI
    print num_II
    
    quality = c_quality + ci_quality + ii_quality
    print quality
    if num_C > 0:
        print c_quality/num_C
    if num_CI > 0:
        print ci_quality/num_CI
    if (num_II > 0):
        print ii_quality/num_II
    results += str(num_matches) + "\t" + str(num_C) + "\t" + str(num_CI) + "\t" + str(num_II) + "\t" + str(quality) + "\t" 
    if num_C > 0:
        results += str(c_quality/num_C) + "\t" 
    else:
        results += "na\t"
    if num_CI > 0:
        results += str(ci_quality/num_CI) + "\t"
    else:
        results += "na\t"
    if (num_II > 0):
        results += str(ii_quality/num_II)
    else:
        results += "na"
    results += "\n"
    pastData.append((num_matches,num_C, num_CI, num_II, quality, 
        c_quality/num_C if num_C else 0, ci_quality/num_CI if num_CI else 0, ii_quality/num_II if num_II else 0))

avgs = np.mean(pastData, axis=0)
stdevs = np.std(pastData, axis=0)

s = ''
for std in stdevs:
    s += str(std)+"\t"
s += "\n"
results = s+"\n\n\n"+results
s = ''
for a in avgs:
    s += str(a)+"\t"
s+="\n"
results = s+results
    
if args.outputFile:
    with open(args.outputFile,'w') as f:
        f.write(results)
else:
    print results
