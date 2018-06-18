from gurobipy import *
import json
import argparse







inputfile = 'data.json'

with open(inputfile, 'r') as f:
    data = json.load(f)


T = data[0]
K = data[1]
matches = data[2]

model = Model("Dual Optimizer")

alpha = {}
beta = {}

for t in range(T+K):
    if sum(matches[t]) > 0:
        alpha[t] = model.addVar(vtype = GRB.CONTINUOUS, lb=0, name='alpha_'+str(t))

for i in range(1,K+1):
    if sum(matches[i+T-1]) > 0:
            beta[i] = model.addVar(vtype=GRB.CONTINUOUS, lb=0, name='beta_'+str(i))
            continue
    for t in range(T+K):
        if matches[t][i] != 0:
            beta[i] = model.addVar(vtype=GRB.CONTINUOUS, lb=0, name='beta_'+str(i))
            break
beta[0] = 0


model.addConstrs((matches[t][i] - alpha[t] - beta[i] - (beta[t-T+1] if t+1-T in beta else 0) <= 0 for t in alpha for i in beta), "something...")

obj = quicksum(alpha[t] for t in alpha) + quicksum(beta[i] for i in beta)
model.setObjective(obj, GRB.MINIMIZE)
model.optimize()
