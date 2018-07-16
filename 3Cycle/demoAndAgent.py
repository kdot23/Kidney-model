#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  3 13:13:30 2018

@author: kelseylieberman
"""
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--data', nargs = "+", help = "List of files with demographic info")
parser.add_argument('--agent', help = "File with agent information")
parser.add_argument('-T', default = 50, type = int)
parser.add_argument('-K', default = 50, type = int)
parser.add_argument('-o)', '--output')
parser.add_argument('--addBeta', help = "Add predicted beta to agent information, betas from agent csv of \
                    onlineMatching of the same population")
args = parser.parse_args()

results = ''
     
with open(args.agent, 'r') as f:
    agent = f.readlines()

if (args.addBeta):
    with open(args.addBeta, 'r') as f:
        beta = f.readlines()

agentList = []
t = 0
for fn in args.data:
    with open(fn, 'rb') as f:
        d = pickle.load(f)
        demo = d[4]
        
    if (args.addBeta):
        betas = {}
        for i in range(t*(args.K + args.T), (t+1)*(args.K + args.T)):
            betaInfo = beta[i].split('\t')
            betas[betaInfo[0]] = betaInfo[6][:-1]
            
    for i in range(t*(args.K+args.T), (t+1)*(args.K+args.T)):
        agentList.append(agent[i].split('\t'))
        
        if (agentList[i][0][0] == 'I'):
            id = int(agentList[i][0][1:]) + args.T-1
        else:
            id = int(agentList[i][0][1:])-1
        if (args.addBeta):           
            results += agent[i][:-1] + "\t" + betas[agentList[i][0]] + "\t" + "\t".join(str(v) for v in demo[id]) + "\n"
        else:            
            #-1 cuts off new line character
            results += agent[i][:-1] + "\t" + "\t".join(str(v) for v in demo[id]) + "\n"
    print t
    t += 1
if args.output:
    with open(args.output, 'w') as f:
        f.write(results)
    
