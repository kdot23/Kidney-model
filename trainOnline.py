import json
from sklearn import linear_model
from sklearn.preprocessing import PolynomialFeatures
import tensorflow as tf
from sklearn.model_selection import cross_val_score
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('--inputFile', default='trainData.json')
parser.add_argument('-i','--input', default='./')
parser.add_argument('--diffPops', action='store_true')
parser.add_argument('-o', '--output')
parser.add_argument('--coef', action='store_true', help = "boolean for coefficients for linear model")
args = parser.parse_args()


varNames = ["bloodTypePatient O", "bloodTypePatient A", "bloodTypePatient B", "bloodTypePatient AB",
            "bloodTypeDonor", "donor_afam", "donor_age", "donor_sex", "donor_cig_use", \
            "rec_sex", "donor_weight", "rec_weight", "donor_bmi"]
varsUsed = [0, 1, 3, 5, 8, 9]

inputFile = args.inputFile

with open(args.input+'/'+inputFile, 'r') as f:
    data = json.load(f)

values = []
labels = []

for d in data:
    demo = d[0]
    values.append([demo[v] for v in varsUsed])
    labels.append(d[1])
results = ''
if args.coef:
    LR = linear_model.LinearRegression()
    LR.fit(values, labels)
    for i in range(len(varsUsed)):
        results += str(LR.coef_[i]) + "\t"
    results += "\n"
"""
#print LR.coef_
#print "intercept: " + str(LR.intercept_)
#print "R^2: " + str(LR.score())
results = ''
"""
if args.coef:
    if args.output:
        with open(args.output, 'w') as f:
            f.write(results)
    else:
        print results
    exit()

if not args.diffPops:
    num_chunks = 10
    chunk_size = len(data) / num_chunks
    for i in range(num_chunks):
        for j in range(1,5):
            test_chunk_val = values[i*chunk_size: (i+1)*chunk_size]
            test_chunk_lab = labels[i*chunk_size: (i+1)*chunk_size]
            train_chunk_val = values[: i*chunk_size] + values[(i+1)*chunk_size:]
            train_chunk_lab = labels[: i*chunk_size] + labels[(i+1)*chunk_size:]
            poly = PolynomialFeatures(degree=j)
            X = poly.fit_transform(train_chunk_val)
            LR = linear_model.LinearRegression()
            LR.fit(X, train_chunk_lab)
            X2 = poly.fit_transform(test_chunk_val)
            results +=  str(LR.score(X2, test_chunk_lab))+"\t"
        results += "\n"
else:
    for fn in os.listdir(args.input):
        if fn == inputFile or fn == '.DS_Store': continue
        data = None
        with open(args.input+'/'+fn, 'r') as f:
            data = json.load(f)
        
        testValues = []
        testLabels = []
        
        for d in data:
            demo = d[0]
            testValues.append([demo[v] for v in varsUsed])
            testLabels.append(d[1])
        for i in range(1,5):
            poly = PolynomialFeatures(degree=i)
            X = poly.fit_transform(values)
            LR = linear_model.LinearRegression()
            LR.fit(X, labels)
            X2 = poly.fit_transform(testValues)
            results += str(LR.score(X2, testLabels)) + "\t"
        results += "\n"


    
if args.output:
    with open(args.output, 'w') as f:
        f.write(results)
else:
    print results


"""
for i in range(2,6):
    poly = PolynomialFeatures(degree=i)
    X = poly.fit_transform(values)
    LR = linear_model.LinearRegression()
    LR.fit(X, labels)
    print "r^2-" + str(i) + ": " + str(LR.score(X,labels))
"""   
