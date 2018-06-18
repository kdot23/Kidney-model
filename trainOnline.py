import json
from sklearn import linear_model

varsUsed = [0, 1, 3, 5, 8, 9]

inputFile = 'trainData.json'

with open(inputFile, 'r') as f:
    data = json.load(f)

values = []
labels = []

for d in data:
    demo = d[0]
    values.append([demo[v] for v in varsUsed])
    labels.append(d[1])

print labels

LR = linear_model.LinearRegression()
LR.fit(values, labels)
print LR.coef_
print LR.intercept_
