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
parser.add_argument('-T', default = 50)
parser.add_argument('-o)', '--output')
args = parser.parse_args()

results = ''
t = 1
for fn in args.data:
    with open(fn, 'rb') as f:
        d = pickle.load(f)
        demo = d[4]
    print t
    t += 1
    agentList = []     
    with open(args.agent, 'r') as f:
        agent = f.readlines()
    for i in range(2*args.T):
        agentList.append(agent[i].split('\t'))
    
    for i in range(len(agentList)):
        if (agentList[i][0][0] == 'I'):
            id = int(agentList[i][0][1:]) + args.T-1
        else:
            id = int(agentList[i][0][1:])-1
        results += agent[i][:-1] + "\t" + "\t".join(str(v) for v in demo[id]) + "\n"

if args.output:
    with open(args.output, 'w') as f:
        f.write(results)
    
