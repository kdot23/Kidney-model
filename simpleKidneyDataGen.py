import random
import json

T = 100
compat = []
#d=0 is more compatible and d=1 is less compatible, r=0 is easy to match, r=1 hard to match
for i in range(T):
    if random.random() < .6:
        d = 0
    else:
        d = 1
    if random.random() < .5:
        r = 1
    else:
        r = 0
    q = random.expovariate(.2)
    compat.append((d,r,q))
print compat
    
incompat = []        
for i in range(T):
    if random.random() < .5:
        d = 0
    else:
        d = 1
    if random.random() < .3:
        r = 1
    else:
        r = 0
    q = random.expovariate(.2)
    incompat.append((d,r,q))
print incompat

Igraph = []
for i in range(T):
    Igraph.append([])
    for j in range(T):
        if i==j:
            Igraph[i].append(0)
            continue
        p = random.random()
        if not incompat[i][0] and not incompat[j][1]:
            Igraph[i].append(int(p < .04))
        if not incompat[i][0] and incompat[j][1]:
            Igraph[i].append(int(p < .03))
        if incompat[i][0] and not incompat[j][1]:
            Igraph[i].append(int(p < .02))
        if incompat[i][0] and incompat[j][1]:
            Igraph[i].append(int(p < .01))
print Igraph

Cgraph = []
for i in range(T):
    Cgraph.append([])
    for j in range(T):
        p = random.random()
        if not compat[i][0] and not incompat[j][1]:
            Cgraph[i].append(int(p < .04))
        if not compat[i][0] and incompat[j][1]:
            Cgraph[i].append(int(p < .03))
        if compat[i][0] and not incompat[j][1]:
            Cgraph[i].append(int(p < .02))
        if compat[i][0] and incompat[j][1]:
            Cgraph[i].append(int(p < .01))
print Cgraph

with open('data.json', 'w') as f:
    f.write(json.dumps((compat, incompat, Igraph, Cgraph)))
