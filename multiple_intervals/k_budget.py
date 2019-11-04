from datetime import datetime, timedelta
import numpy as np
import copy
import sys
import pickle
import k_utils as utils

def runKBudget(timestamps, k, maxiters = 10, m = {}):
    """
    Implements k-Budget algorithm

    Parameters
    ----------
    timestamps : list of tuples
        Sorted list of interactions [(t1, n1, n2), (t2, n3, n4),...], where t is timestamp, n1 and n2 are interactiong nodes.
        Nodes in the interactions are sorted lexicographically.
    k : int 
        number of activity intervals
    maxiters : int
        maximum number of interactions in binary search and gap points search
    m : dict
        dict of gap points {n: [t1, t2, ..]} where n is a node, list t1, t2 .. is a list of its k-1 inactive (gap) points

    Returns
    -------
    tuple of dicts
        (dict1, dict2): dict1 of a shape {n1: {0: t10, 1: t11, 2: t12}, n2: {0: t20, 1: t21, 2: t22}, ..}, tij is a starting point of jth activity interval of node ni.
                        dict2 of a shape {n1: {0: t10, 1: t11, 2: t12}, n2: {0: t20, 1: t21, 2: t22}, ..}, tij is an ending point of jth activity interval of node ni.

    """
    nodeEdgeIndex = utils.getIndex(timestamps, 'budget')

    if not m:
        m = getGapPoints(nodeEdgeIndex, timestamps[0][0],timestamps[-1][0], k)
    Xstart,  Xend = runKBudget_BS(timestamps, nodeEdgeIndex, k, maxiters = maxiters, m = m)
    for i in xrange(maxiters):
        m = getNewGapPoints(nodeEdgeIndex, timestamps[0][0],timestamps[-1][0], Xstart,  Xend)
        Xstart,  Xend = runKBudget_BS(timestamps, nodeEdgeIndex, k, maxiters = maxiters, m = m)
        
    return Xstart,  Xend

def getGapPoints(nodeEdgeIndex, left_border,right_border, k):
    m = {}    
    for n, ts in nodeEdgeIndex.iteritems():
        s = [left_border] + sorted(ts.keys()) + [right_border]
        
        gaps =  np.argsort([(s[i+1] - s[i]) for i in xrange(len(s)-1)])[::-1]        
        m[n] = sorted([(s[i]+s[i+1])/2.0 for i in gaps[:min((k-1),len(gaps))]])
    return m
    
def getTrueGapPoints(nodeEdgeIndex, active_truths):
    m = {}  
    for n, ints in active_truths.iteritems():
        m[n] = []
        a = ints[0][1]
        for i in xrange(1,len(ints)):
            b = ints[i][0]            
            m[n].append((a + b)/2.0)
            a = ints[i][1]
    return m
    
def getNewGapPoints(nodeEdgeIndex, left_border, right_border, Xstart,  Xend):
    m = {}  
    for n in Xstart.keys():
        m[n] = []
        a = Xend[n][0]
        for i in xrange(1,len(Xstart[n])):
            b = Xstart[n][i]
            if a == -np.Inf:
                a = float(left_border)
            if b == np.Inf:
                b = float(right_border)
            m[n].append((a + b)/2.0)
            a = Xend[n][i]
    return m

def setbudget(nodeEdgeIndex, kfactor):
    b = {}    
    for n, ts in nodeEdgeIndex.iteritems():
        b[n] = (max(ts.keys()) - min(ts.keys()))*kfactor
    return b
    
def DFS(v, t, nodeEdgeIndex, O, Orev, b, S):
    S.append((v,t))
    nodeEdgeIndex[v][t]['visited'] = True    
    
    i = Orev[v][t]
    while O[v][i] and O[v][i][0][0] < (t-b[v][i]):
        (s,u) = O[v][i].pop(0)        
        if not nodeEdgeIndex[u][s]['visited']:
            DFS(u, s, nodeEdgeIndex, O, Orev, b, S)
        
    while O[v][i] and O[v][i][-1][0] > (t+b[v][i]):
        (s,u) = O[v][i].pop(-1)
        if not nodeEdgeIndex[u][s]['visited']:
            DFS(u, s, nodeEdgeIndex, O, Orev, b, S)
    return 
    
def DFS_iterative(v, t, nodeEdgeIndex, O, Orev, b, S):
    
    localS = []
    localS.append((v,t))
    while localS:
        (v,t) = localS.pop()
        if not nodeEdgeIndex[v][t]['visited']:
            nodeEdgeIndex[v][t]['visited'] = True
            S.append((v,t))
            i = Orev[v][t]
            
            tmp = []
            if O[v][i]:
                counter = 0
                while O[v][i][counter][0] < (t-b[v][i]):
                    tmp.append((O[v][i][counter][1], O[v][i][counter][0]))
                    counter += 1
                counter = -1
                while O[v][i][counter][0] > (t+b[v][i]):
                    tmp.append((O[v][i][counter][1], O[v][i][counter][0]))
                    counter -= 1               
            
            for pair in tmp[::-1]:              
                localS.append(pair)
         
                
    return 
    
def reverseDFS_iterative(v, t, nodeEdgeIndex, Q, Orev, b ,R):
    localS = []
    localS.append((v,t))
    
    while localS:
        (v,t) = localS.pop()
        if not nodeEdgeIndex[v][t]['visited']:
            nodeEdgeIndex[v][t]['visited'] = True
            R.append((v,t))    
            for (n1,n2) in nodeEdgeIndex[v][t]['edges']:
                u = n2 if n1 == v else n1
                i = Orev[u][t]
                
                tmp = []
                if Q[u][i]:
                    counter = 0
                    while Q[u][i][counter] < (t-b[u][i]):
                        tmp.append((u, Q[u][i][counter]))
                        counter += 1
                    counter = -1
                    while Q[u][i][counter] > (t+b[u][i]):
                        tmp.append((u, Q[u][i][counter]))
                        counter -= 1
                    
                for pair in tmp[::-1]:
                    localS.append(pair)
    return 
    
def reverseDFS(v, t, nodeEdgeIndex, Q, Orev, b ,R):
    R.append((v,t))
    nodeEdgeIndex[v][t]['visited'] = True
    for (n1,n2) in nodeEdgeIndex[v][t]['edges']:
        u = n2 if n1 == v else n1
        i = Orev[u][t]
        
        while Q[u][i] and Q[u][i][0] < (t-b[u][i]):
            s = Q[u][i].pop(0)
            if not nodeEdgeIndex[u][s]['visited']:
                reverseDFS(u, s, nodeEdgeIndex, Q, Orev, b, R)
                
        while Q[u][i] and Q[u][i][-1] > (t+b[u][i]):
            s = Q[u][i].pop(-1)
            if not nodeEdgeIndex[u][s]['visited']:
                reverseDFS(u, s, nodeEdgeIndex, Q, Orev, b, R)

    return 
                   
def clearVisited(nodeEdgeIndex):
    for n, ts in nodeEdgeIndex.iteritems():       
        for t in ts:
            ts[t]['visited'] = False
    return
    
def kbudgetAlgorithm(timestamps, nodeEdgeIndex, b, m):
    clearVisited(nodeEdgeIndex)
    O, Orev = getO(timestamps, m)
    Smain = []
    for n, ts in nodeEdgeIndex.iteritems():       
        for t in ts:
            if not ts[t]['visited']:
                DFS_iterative(n, t, nodeEdgeIndex, O, Orev,  b, Smain)
    clearVisited(nodeEdgeIndex)
    Q = getQ(timestamps, m)
    R = []
    for (node, tm) in Smain[::-1]:
        if not nodeEdgeIndex[node][tm]['visited']:
            reverseDFS_iterative(node, tm, nodeEdgeIndex, Q, Orev, b, R)   
    
    p1 = {n: { ind: -np.Inf for ind in xrange(len(mpoints)+1)} for n, mpoints in m.iteritems()}
    i2 = {n: { ind: -np.Inf for ind in xrange(len(mpoints)+1)} for n, mpoints in m.iteritems()}
    p2 = {n: { ind: np.Inf for ind in xrange(len(mpoints)+1)} for n, mpoints in m.iteritems()}
    i1 = {n: { ind: np.Inf for ind in xrange(len(mpoints)+1)} for n, mpoints in m.iteritems()}
  
    for (v,t) in R[::-1]:
        i = Orev[v][t]
        if p1[v][i] <= t and p2[v][i] >= t:
            p1[v][i] = max(p1[v][i], t - b[v][i])
            p2[v][i] = min(p2[v][i], t + b[v][i])
            i1[v][i] = min(i1[v][i], t)
            i2[v][i] = max(i2[v][i], t)
            
    Xstart = i1
    Xend = i2
    return Xstart, Xend

 
def getO(timestamps, m):   
    O = {n: { ind:[] for ind in xrange(len(mpoints)+1)} for n, mpoints in m.iteritems()}
    ORev = {n: {} for n in m}
    counter = {n: 0 for n in m}
    
    for (t, n1, n2) in timestamps:
    
        if counter[n1] < len(m[n1]) and t > m[n1][counter[n1]]:
            counter[n1] += 1            
        O[n1][counter[n1]].append((t, n2))        
        ORev[n1][t] = counter[n1]    
        
        if counter[n2] < len(m[n2]) and t > m[n2][counter[n2]]:
            counter[n2] += 1
        O[n2][counter[n2]].append((t,n1))
        ORev[n2][t] = counter[n2]
    return O, ORev
        
def getQ(timestamps, m):
    O = {n: { ind: [] for ind in xrange(len(mpoints)+1)} for n, mpoints in m.iteritems()}
    counter = {n: 0 for n in m}
    
    for (t, n1, n2) in timestamps:
        if counter[n1] < len(m[n1]) and t > m[n1][counter[n1]]:
            counter[n1] += 1
        O[n1][counter[n1]].append(t)        
                
        if counter[n2] < len(m[n2]) and t > m[n2][counter[n2]]:
            counter[n2] += 1
        O[n2][counter[n2]].append(t)
        O = {n: {i: sorted(set(l)) for i, l in ind.iteritems()} for n, ind in O.iteritems()}
        #O = {n: {i: sorted(set(l), reverse=True) for i, l in ind.iteritems()} for n, ind in O.iteritems()}
    return O
    
def runKBudget_fixed_budget(timestamps, nodeEdgeIndex, b = {}, m = {}, k = 10, kfactor = 1.0):
    if not b:
        b = {n: [(timestamps[-1][0]-timestamps[0][0])]*k for n in nodeEdgeIndex}
    if not m:
        m = getGapPoints(nodeEdgeIndex, timestamps[0][0],timestamps[-1][0], k)
    Xstart, Xend = kbudgetAlgorithm(timestamps, nodeEdgeIndex, b, m)
   
    return Xstart, Xend
    

    
    
def runKBudget_BS(timestamps, nodeEdgeIndex, kintervals, klow = -1, kup = -1, maxiters = 10, m = {}):
    
    maxb = timestamps[-1][0]-timestamps[0][0]
    if klow ==  -1:
        klow, kup = 1./maxb, 1.0 
    b = {n: [maxb]*kintervals for n in nodeEdgeIndex}
    if not m:
        m = getGapPoints(nodeEdgeIndex, timestamps[0][0],timestamps[-1][0], kintervals)
    LXstart, LXend = kbudgetAlgorithm(timestamps, nodeEdgeIndex, b, m)
    
    alpha = 1.0/maxb
    itercount = 0
    while (kup - klow) > alpha and itercount < maxiters:
        #print 'iteration', itercount, kup - klow, alpha
        itercount += 1
        k = (kup + klow)/2
        b = {n: [maxb*k]*kintervals for n in nodeEdgeIndex}
        Xstart, Xend = kbudgetAlgorithm(timestamps, nodeEdgeIndex, b, m)
        uncovered =  utils.checkCoverage(Xstart, Xend, timestamps)
        if uncovered:
            klow = k
        else:
            kup = k
            LXstart, LXend = Xstart, Xend
            
    return LXstart, LXend
