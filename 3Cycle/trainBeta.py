import argparse
import json
import pickle
from gurobipy import *
import numpy as np
import random
from sklearn import linear_model
from sklearn import ensemble
from sklearn.preprocessing import PolynomialFeatures
import os
import matplotlib.pyplot as plt
plt.style.use('ggplot')
import pdb

cl1 = ["#006D2C", "#31A354", "#74C476"]

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


def calcBetaLP(T, K, matches, used_incompat):
    estimator = Model("estimate beta vals")
    alpha = {}
    beta = {}
    for i in range(1,K+1):
        if i in used_incompat: continue
        if any(k[0] == i+T for k in matches if k[1] not in used_incompat and k[2] not in used_incompat):
            beta[i] = estimator.addVar(vtype=GRB.CONTINUOUS, lb=0, name='beta_'+str(i))
        elif any(k[1] == i for k in matches if k[0] > T and k[0]-T not in used_incompat and k[2] not in used_incompat):
            beta[i] = estimator.addVar(vtype=GRB.CONTINUOUS, lb=0, name='beta_'+str(i))
        elif any(k[2] == i for k in matches if k[0] > T and k[0]-T not in used_incompat and k[1] not in used_incompat):
            beta[i] = estimator.addVar(vtype=GRB.CONTINUOUS, lb=0, name='beta_'+str(i))
    if args.quality:
        estimator.addConstrs((matches[t,i,j] - beta[i] -beta[j] - (beta[t-T] if t-T in beta else 0) <= 0 for t in range(T+1,T+K+1) if t-T in beta for i in beta \
                for j in beta if (t,i,j) in matches), 'something')
    else:
        estimator.addConstrs((COUNT((t,i,j)) - beta[i] - beta[j] - (beta[t-T] if t-T in beta else 0) <= 0 for t in range(T+1,T+K+1) if t-T in beta for i in beta \
                for j in beta if (t,i,j) in matches), 'something')
    obj = quicksum(beta[i] for i in beta)
    estimator.setObjective(obj, GRB.MINIMIZE)
    estimator.optimize()
    newBeta = {i:beta[i].X for i in beta}
    for i in range(1, K+1):
        if i in used_incompat:continue
        if i not in newBeta:
            newBeta[i] = 0
    return newBeta

parser = argparse.ArgumentParser()
parser.add_argument('--trainFiles', nargs = "+", help = "List of files to train on")
parser.add_argument('--testFiles', nargs = "+", help = "List of files to test")
parser.add_argument('-o', '--output', help = 'csv file to output count and quality to')
parser.add_argument('--agents', help='output the quality of each agent to this file (.csv)')
parser.add_argument("-v", "--useVars", nargs = "+", type = int, default=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18], \
                    help = "List of variables")
parser.add_argument('--quality', action='store_true', help='Flag should be present if optimization is being done for quality')
parser.add_argument('--forestRegression', nargs='?', const=10, type=int, \
        help='Flag should be present if forest regression is to be used instead of Linear, optional argument of number of trees')
parser.add_argument('--graph', help='stem of output file for graph')
parser.add_argument('-d', '--degree', default=1, type=int, help='type of polynomial to use while training')
parser.add_argument('--lpEstimator', action='store_true', help='should be present if dual on incompatibles only should be used to \
        estimate betas')
parser.add_argument('--lpRepeat', action='store_true', help='should be present if lp repeat method is used to estimate betas')
parser.add_argument('--graph_state', action='store_true', help='flag should be present if lpEstimator values are going to be used for training')
args = parser.parse_args()


varsUsed = args.useVars
data = []
if args.trainFiles:
    for fn in args.trainFiles:
        with open(fn, 'r') as f:
            data += json.load(f)
    random.shuffle(data)

    values = []
    labels = []

    for d in data:
        demo = d[0]
        if args.graph_state:
            values.append([demo[v] for v in varsUsed] + [d[2]])
        else:
            values.append([demo[v] for v in varsUsed])
        labels.append(d[1])

    # poly = PolynomialFeatures(degree=args.degree)
    X = np.array(values)
    if not args.forestRegression:
        LR = linear_model.Ridge()
    else:
        LR = ensemble.RandomForestRegressor(n_estimators=args.forestRegression)

    LR.fit(X, labels)

dataIndex = 0
data = []
for fn in args.testFiles:
    with open(fn, 'r') as f:
        data += json.load(f)

    values = []
    labels = []

    for d in data:
        demo = d[0]
        if args.graph_state:
            values.append([demo[v] for v in varsUsed] + [d[2]])
        else:
            values.append([demo[v] for v in varsUsed])
        labels.append(d[1])

    # poly = PolynomialFeatures(degree=args.degree)
    X2 = np.array(values)
    predicted_labels = LR.predict(X2)

importances = LR.feature_importances_
print importances

error = predicted_labels - np.array(labels)
final_error = np.sum(error * error)
print final_error

std = np.std([tree.feature_importances_ for tree in LR.estimators_],
             axis=0)
indices = np.argsort(importances)[::-1]

# Print the feature ranking
print("Feature ranking:")

for f in range(X.shape[1]):
    print("%d. feature %d (%f)" % (f + 1, indices[f], importances[indices[f]]))

# Plot the feature importances of the forest
feature_name = ["donor O", "donor A", "donor B", "donor AB",
                "rec O", "rec A", "rec B", "rec AB",
                "donor AA", "donor age",
                "donor sex", "donor cig use", "rec sex",
                "donor weight", "rec weight", "donor bmi", "donor egfr",
                "donor sbp", "donor PRA", "beta"]



plt.figure()
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.size"] = 20

# plt.title("Feature importances")
plt.bar(range(X.shape[1]), importances[indices], yerr=std[indices], align="center")
plt.xticks(range(X.shape[1]), [feature_name[i] for i in indices],rotation=90)
plt.xlim([-1, X.shape[1]])

ax = plt.gca()
ax.spines['bottom'].set_color('black')
ax.spines['left'].set_color('black')
ax.spines['top'].set_color('black')
ax.spines['right'].set_color('black')
ax.tick_params(axis='x', colors='black')
ax.tick_params(axis='y', colors='black')
ax.set_axis_bgcolor('white')

plt.tight_layout()

plt.savefig("feature_importance.pdf")

plt.figure()
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.size"] = 20
plt.xlabel(r"True $\beta$", color="black")
plt.ylabel(r"Predicted $\beta$", color="black")
plt.scatter(labels, predicted_labels)
plt.plot(labels, labels, color='black')
# plt.legend()

ax = plt.gca()
ax.spines['bottom'].set_color('black')
ax.spines['left'].set_color('black')
ax.spines['top'].set_color('black')
ax.spines['right'].set_color('black')
ax.tick_params(axis='x', colors='black')
ax.tick_params(axis='y', colors='black')
ax.set_axis_bgcolor('white')

plt.tight_layout()
plt.savefig("true_predicted_scatter.pdf")
plt.show()
    # with open(fn, 'r') as f:
    #     data = pickle.load(f)
    #
    # K = data[0]
    # T = data[1]
    # matches = data[3]
    # demo = data[4]
    # directed_matches = data[6]
    # used_incompat = set()
    #
    # if args.graph_state:
    #     lpbeta = calcBetaLP(T, K, matches, used_incompat)
    #     testValues = [[demo[i][v] for v in varsUsed] + [lpbeta[i - T + 1]] for i in range(T, T + K)]
    # else:
    #     testValues = [[demo[i][v] for v in varsUsed] for i in range(T, T + K)]
    # X2 = poly.fit_transform(testValues)
    #
    # if not (args.lpEstimator or args.lpRepeat):
    #     X2 = poly.fit_transform(testValues)
    #     betaList = LR.predict(X2)
    #     beta = {i + 1: betaList[i] for i in range(len(betaList))}
    # else:
    #     beta = calcBetaLP(T, K, matches, used_incompat)


