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
from SaidmanCompatibleGenerator import BJCSensitivityPool, DistributionGenerator
import functions
from DistributionGenerator import *
import util



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
parser.add_argument('--graph_state', action='store_true', help='flag should be present if lpEstimator values are going to be used for training')
parser.add_argument('--fwd_proj', nargs='?', const=10, help='flag should be present if the simulator uses forward projections to estimate beta values')
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

def calcBetasProj(T, K, matches):
    estimator = Model('estimate beta values')
    alpha = {}
    beta = {}
    for i in range(1, T+1):
        if any(v[0] == i for v in matches):
            alpha[i] = estimator.addVar(vtype = GRB.CONTINUOUS, lb = 0, name = 'alpha_'+str(i))
    for i in range(1, K+1):
        if any(v[0] == i+T or v[1] == i for v in matches):
            beta[i] = estimator.addVar(vtype = GRB.CONTINUOUS, lb = 0, name = 'beta+'+str(i))
    beta[0] = 0
    if args.quality:
        estimator.addConstrs((matches[t,i] - (alpha[t] if t in alpha else 0) -  beta[i] -  (beta[t-T] if t-T in beta else 0) <= 0 \
                for t in range(1,T+K+1) if t-T in beta or t in alpha for i in beta if (t,i) in matches),  'something...')
    else:
        estimator.addConstrs((COUNT((t,i)) - (alpha[t] if t in alpha else 0) - beta[i] - (beta[t-T] if t-T in beta else 0) <= 0 \
                for t in range(1, T+K+1) if t-T in beta or t in alpha for i in beta if (t,i) in matches), 'something...')
    obj = quicksum(beta[i] for i in beta) + quicksum(alpha[t] for t in alpha)
    estimator.setObjective(obj, GRB.MINIMIZE)
    estimator.optimize()
    newBeta =  {i:(beta[i].X if i in beta else 0) for i in range(1,K+1)}
    for i in range(0, K+1):
        if i in used_incompat:
            continue
        if i not in newBeta:
            newBeta[i] = 0
    return newBeta

def calcBetasLP(T, K, matches, used_incompat):
    estimator = Model('estimate beta values')
    alpha = {}
    beta = {}
    for i in range(1,K+1):
        if i in used_incompat: continue
        if any(k[0] == i+T for k in matches if k[1] not in used_incompat) or any(k[1] == i for k in matches if k[0] > T \
                and k[0]-T not in used_incompat):
            beta[i] = estimator.addVar(vtype=GRB.CONTINUOUS, lb=0, \
                    name='beta_'+str(i))
    if args.quality:
        estimator.addConstrs((matches[t,i] -  beta[i] -  (beta[t-T] if t-T in beta else 0) <= 0 \
                for t in range(T+1,T+K+1) if t-T in beta for i in beta if (t,i) in matches),  'something...')
    else:
        estimator.addConstrs((COUNT((t,i)) - beta[i] - (beta[t-T] if t-T in beta else 0) <= 0 \
                for t in range(T+1, T+K+1) if t-T in beta for i in beta if (t,i) in matches), 'something...')
    obj = quicksum(beta[i] for i in beta)
    estimator.setObjective(obj, GRB.MINIMIZE)
    estimator.optimize()
    newBeta =  {i:beta[i].X for i in beta}
    for i in range(1, K+1):
        if i in used_incompat:
            continue
        if i not in newBeta:
            newBeta[i] = 0
    return newBeta


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
    if args.graph_state:
        values.append([demo[v] for v in varsUsed]+[d[2]])
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

    K = data[0]
    T = data[1]
    matches = data[2]
    demo = data[4]
    directed_matches = data[6]
    IMisMatches = data[5]
    used_incompat = set()
    
    testValues = [[demo[i][v] for v in varsUsed] for i in range(T,T+K)]
    if args.graph_state:
        lpBeta = calcBetasLP(T,K,matches,used_incompat)
        for i in range(len(testValues)):
            testValues[i].append(lpBeta[i+1])
            
    X2 = poly.fit_transform(testValues)
    if not (args.lpEstimator or args.lpRepeat or args.fwd_proj):
        betaList = LR.predict(X2)
        beta = {i+1:betaList[i] for i in range(len(betaList))}
    elif args.fwd_proj:
        betas = []
        gen = DistributionGenerator()
        Imatches = {v:matches[v] for v in matches if v[0] > T}
        def getBloodTypeDonor(i):
            if demo[i][0]: return 0
            elif demo[i][1]: return 1
            elif demo[i][2]: return 2
            else: return 3
        def getBloodTypeRecip(i):
            if demo[i][4]: return 0
            elif demo[i][5]: return 1
            elif demo[i][6]: return 2
            else: return 3
        for _ in range(args.fwd_proj):
            proj_matches = dict(Imatches)
            pool = BJCSensitivityPool(T, 0)
            Tbeta = []
            misMatches = {}
            positiveCrossMatches = {}
            for i in range(T+K):
                for j in range(T+K):
                    if j < T and i < T: continue
                    misMatches[i+1, j+1] = (gen.gen_donor_rec_HLA_B_mis(0), gen.gen_donor_rec_HLA_DR_mis(0))
            for t in pool.compatiblePairs:
                for i in range(T+1,T+K+1):
                    positiveCrossMatches[t, i] = t.saidman.isPositiveCrossmatch(demo[i-1][18])
                    positiveCrossMatches[i, t] = t.saidman.isPositiveCrossmatch(t.pr_PraIncompatiblity)
            Cset = set(pool.compatiblePairs)
            def compatible(donor, recipient):
                if donor in Cset:
                    return functions.are_blood_compatible(donor.bloodTypeDonor, getBloodTypeRecip(recipient-1)) \
                              and (not positiveCrossMatches[donor, recipient])
                elif recipient in Cset:
                    return functions.are_blood_compatible(getBloodTypeDonor(donor-1), recipient.bloodTypePatient) \
                            and (not positiveCrossMatches[donor, recipient])
                else:
                    return functions.are_blood_compatible(getBloodTypeDonor(donor-1), getBloodTypeRecip(recipient-1)) \
                            and (not positiveCrossMatches[donor, recipient])
            def getLKDPI(donor, recipient, HLA_B_mis, HLA_DR_mis ):
                if donor in Cset:
                    bloodTypeDonor = donor.bloodTypeDonor
                    donor_weight = donor.donor_weight
                    donor_age = donor.donor_age
                    donor_afam = donor.donor_afam
                    donor_bmi = donor.donor_bmi
                    donor_cig_use = donor.donor_cig_use
                    donor_sex = donor.donor_sex
                    donor_sbp = donor.donor_sbp
                    donor_egfr = donor.donor_egfr
                else:
                    bloodTypeDonor = getBloodTypeDonor(donor-1)
                    donor_weight = demo[donor-1][13]
                    donor_age = demo[donor-1][9]
                    donor_afam = demo[donor-1][8]
                    donor_bmi = demo[donor-1][15]
                    donor_cig_use = demo[donor-1][11]
                    donor_sex = demo[donor-1][10]
                    donor_sbp = demo[donor-1][17]
                    donor_egfr = demo[donor-1][16]
                if recipient in Cset:
                    bloodTypePatient = recipient.bloodTypePatient
                    rec_weight = recipient.rec_weight
                    rec_sex = recipient.rec_sex
                else:
                    bloodTypePatient = getBloodTypeRecip(recipient-1)
                    rec_weight = demo[recipient-1][14]
                    rec_sex = demo[recipient-1][12]
                donor_rec_abo_comp = functions.are_blood_compatible(bloodTypeDonor, bloodTypePatient)
                weight_ratio = util.calculate_weight_ratio(donor_weight, rec_weight)
                return util.calculate_lkdpi(donor_age, donor_afam, donor_bmi, donor_cig_use, \
                                                  donor_sex, rec_sex, donor_sbp, donor_rec_abo_comp, \
                                                  0, donor_egfr, HLA_B_mis, HLA_DR_mis, weight_ratio)
            for t in range(1,T+1):
                for i in range(1, K+1):
                    if compatible(pool.compatiblePairs[t-1], i+T) and compatible(i+T, pool.compatiblePairs[t-1]):
                        proj_matches[t,i] = util.calculate_survival(getLKDPI(pool.compatiblePairs[t-1], i, misMatches[t,i+T][0], misMatches[t,i+T][1])) \
                                + util.calculate_survival(getLKDPI(i, pool.compatiblePairs[t-1], misMatches[i+T,t][0], misMatches[i+T,t][1]))
            betaProj = calcBetasProj(T, K, proj_matches)
            betas.append(betaProj)
        beta = {i:sum(b[i] for b in betas)/len(betas) for i in range(K+1)}

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
    
if args.output:
    with open(args.output, 'w') as f:
        f.write(results)
else:
    print results
    
if args.agents:
    with open(args.agents, 'w') as f:
        f.write(agentInfo)   
