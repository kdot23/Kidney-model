import json
from sklearn import linear_model
from sklearn.preprocessing import PolynomialFeatures
from sklearn.model_selection import cross_val_score
import argparse
import os
import random
import pickle

parser = argparse.ArgumentParser()
parser.add_argument('--inputFiles', default=['trainData.json'], nargs='+', help='list of files to use as input, \
if diffpops is present these will be used as train files and all other files in input directory will be used as test files, \
otherWise n fold cross validation is performed')
parser.add_argument('-i','--input', default='./', help='input directory to be used')
parser.add_argument('--diffPops', action='store_true', help='flag for whether or not to use cross validation, see inputFiles flag')
parser.add_argument('-d', '--degree', default=[1], type=int, nargs='+', help='types of polynomials to use while training')
parser.add_argument('-v', '--useVars', default=[0,1,3,5,8,9], type=int, nargs='+', help='variables to train off of' )
parser.add_argument('-o', '--output', help='output File')
parser.add_argument('-n', '--numChunks', default=10, type=int, 'n for n-fold cross validation')
parser.add_argument('--coef', action='store_true', help = "boolean for coefficients for linear model")
parser.add_argument('--saveModel', help='optional flag to save model, argument is stem of file(e.g. if file is given, then model of degree 3 \
        will be saved to file3.dat)')
args = parser.parse_args()


varNames = ["bloodTypePatient", "bloodTypeDonor", "donor_afam", "donor_age", "donor_sex", "donor_cig_use", \
            "rec_sex", "donor_weight", "rec_weight", "donor_bmi"]

data = []
for fn in args.inputFiles:
    with open(args.input+'/'+fn, 'r') as f:
        data += json.load(f)

random.shuffle(data)

values = []
labels = []

for d in data:
    demo = d[0]
    values.append([demo[v] for v in args.useVars])
    labels.append(d[1])
results = ''



if args.coef:
    LR = linear_model.LinearRegression()
    LR.fit(values, labels)
    for i in range(len(args.useVars)):
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
    chunk_size = len(data) / args.numChunks
    for i in range(args.numChunks):
        for j in args.degree:
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
        if fn in args.inputFiles or fn == '.DS_Store': continue
        data = None
        with open(args.input+'/'+fn, 'r') as f:
            data = json.load(f)
        
        testValues = []
        testLabels = []
        
        for d in data:
            demo = d[0]
            testValues.append([demo[v] for v in args.useVars])
            testLabels.append(d[1])
        for i in args.degree:
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

if args.saveModel:
    for i in args.degree:
        poly = PolynomialFeatures(degree=i)
        X = poly.fit_transform(values)
        LR = linear_model.LinearRegression()
        LR.fit(X,labels)
        with open(args.saveModel+str(i)+'.dat', 'wb') as f:
            pickle.dump(LR, f)


"""
for i in range(2,6):
    poly = PolynomialFeatures(degree=i)
    X = poly.fit_transform(values)
    LR = linear_model.LinearRegression()
    LR.fit(X, labels)
    print "r^2-" + str(i) + ": " + str(LR.score(X,labels))
"""   
