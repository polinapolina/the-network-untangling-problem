from datetime import datetime, timedelta
import numpy as np
import copy
import sys
import pickle
import utils

sys.setrecursionlimit(10000000)

def runBudget(timestamps, maxiter = 10):
    """
    Implements Budget algorithm

    Parameters
    ----------
    timestamps : list of tuples
        Sorted list of interactions [(t1, n1, n2), (t2, n3, n4),...], where t is timestamp, n1 and n2 are interactiong nodes.
        Nodes in the interactions are sorted lexicographically.
    maxiter : int
        maximum number of interactions in binary search

    Returns
    -------
    tuple of dicts
        (dict1, dict2): dict1 of a shape {n1: t1, n2: t2, ..} contains starting point t1 of activity interval of node n1, for every node,
                        dict2 of a shape {n1: t1, n2: t2, ..} contains ending point t1 of activity interval of node n1, for every node.

    """
    
    nodeEdgeIndex = utils.getIndex(timestamps, 'budget')
    maxb = timestamps[-1][0]-timestamps[0][0]
    klow, kup = 1./maxb, 1.0 
    b = {n: maxb for n in nodeEdgeIndex}
    LXstart, LXend = budgetAlgorithm(timestamps, nodeEdgeIndex, b)
    c = 1
    
    alpha = 1.0/maxb
    while (kup - klow) > alpha and c < maxiter:
        k = (kup + klow)/2
        b = {n: maxb*k for n in nodeEdgeIndex}
        Xstart, Xend = budgetAlgorithm(timestamps, nodeEdgeIndex, b)
        c += 1
        uncovered =  utils.checkCoverage(Xstart, Xend, timestamps)
        if uncovered:
            klow = k
        else:
            kup = k
            LXstart, LXend = Xstart, Xend
    return LXstart, LXend






def setbudget(nodeEdgeIndex, K):
    b = {}    
    for n, ts in nodeEdgeIndex.iteritems():
        b[n] = (max(ts.keys()) - min(ts.keys()))*K
    return b
    
def DFS(v, t, nodeEdgeIndex, O, b,S):
    S.append((v,t))
    nodeEdgeIndex[v][t]['visited'] = True
    
    while O[v] and O[v][0][0] < (t-b[v]):        
        (s,u) = O[v].pop(0)        
        if not nodeEdgeIndex[u][s]['visited']:
            DFS(u, s, nodeEdgeIndex,O, b, S)
        
    while O[v] and O[v][-1][0] > (t+b[v]):
        (s,u) = O[v].pop(-1)
        if not nodeEdgeIndex[u][s]['visited']:
            DFS(u, s, nodeEdgeIndex, O, b, S)
    return 


def DFS_iterative(v, t, nodeEdgeIndex, O, b, S):
    
    localS = []
    localS.append((v,t))
    while localS:
        (v,t) = localS.pop()
        if not nodeEdgeIndex[v][t]['visited']:
            nodeEdgeIndex[v][t]['visited'] = True
            S.append((v,t))
            
            tmp = []
            if O[v]:
                counter = 0
                while O[v][counter][0] < (t-b[v]):
                    tmp.append((O[v][counter][1], O[v][counter][0]))
                    counter += 1
                counter = -1
                while O[v][counter][0] > (t+b[v]):
                    tmp.append((O[v][counter][1], O[v][counter][0]))
                    counter -= 1               
            
            for pair in tmp[::-1]:              
                localS.append(pair)
                
    return 
    
    
def reverseDFS(v, t, nodeEdgeIndex, Q, b ,R):
    R.append((v,t))
    nodeEdgeIndex[v][t]['visited'] = True
    sorted_t = sorted(nodeEdgeIndex[v].keys())
    for (n1,n2) in nodeEdgeIndex[v][t]['edges']:
        u = n2 if n1 == v else n1
        while Q[u] and Q[u][0] < (t-b[u]):
            s = Q[u].pop(0)        
            if not nodeEdgeIndex[u][s]['visited']:
                reverseDFS(u, s, nodeEdgeIndex, Q, b, R)
                
        while Q[u] and Q[u][-1] > (t+b[u]):
            s = Q[u].pop(-1)        
            if not nodeEdgeIndex[u][s]['visited']:
                reverseDFS(u, s, nodeEdgeIndex, Q, b, R)

    return 

def reverseDFS_iterative(v, t, nodeEdgeIndex, Q, b ,R):
    localS = []
    localS.append((v,t))
    
    while localS:
        (v,t) = localS.pop()
        if not nodeEdgeIndex[v][t]['visited']:
            nodeEdgeIndex[v][t]['visited'] = True
            R.append((v,t))    
            for (n1,n2) in nodeEdgeIndex[v][t]['edges']:
                u = n2 if n1 == v else n1
                
                tmp = []
                if Q[u]:
                    counter = 0
                    while Q[u][counter] < (t-b[u]):
                        tmp.append((u, Q[u][counter]))
                        counter += 1
                    counter = -1
                    while Q[u][counter] > (t+b[u]):
                        tmp.append((u, Q[u][counter]))
                        counter -= 1
                    
                for pair in tmp[::-1]:
                    localS.append(pair)
    return 
                   
def clearVisited(nodeEdgeIndex):
    for n, ts in nodeEdgeIndex.iteritems():       
        for t in ts:
            ts[t]['visited'] = False
    return
    
def budgetAlgorithm(timestamps, nodeEdgeIndex, b):
    clearVisited(nodeEdgeIndex)
    O = getO(timestamps)
    Smain = []
    for n, ts in nodeEdgeIndex.iteritems():       
        for t in ts:
            if not ts[t]['visited']:
                #DFS(n, t, nodeEdgeIndex, O,  b, Smain)
                DFS_iterative(n, t, nodeEdgeIndex, O,  b, Smain)
    clearVisited(nodeEdgeIndex)
    Q = getQ(timestamps)
    R = []
    for (node, tm) in Smain[::-1]:
        if not nodeEdgeIndex[node][tm]['visited']:
            #reverseDFS(node, tm, nodeEdgeIndex, Q, b ,R)
            reverseDFS_iterative(node, tm, nodeEdgeIndex, Q, b ,R)
                
    p1 = {n: -np.Inf for n in nodeEdgeIndex}
    i2 = {n: -np.Inf for n in nodeEdgeIndex}
    p2 = {n: np.Inf for n in nodeEdgeIndex}
    i1 = {n: np.Inf for n in nodeEdgeIndex}   
    
    for (v,t) in R[::-1]:
        if p1[v] <= t and p2[v] >= t:
            p1[v] = max(p1[v], t - b[v])
            p2[v] = min(p2[v], t + b[v])
            i1[v] = min(i1[v], t)
            i2[v] = max(i2[v], t)
            
    Xstart = i1
    Xend = i2
    return Xstart, Xend

 
def getO(timestamps):
    O = {}
    for (t,n1,n2) in timestamps:
        if n1 not in O:
            O[n1] = []
        O[n1].append((t,n2))
        if n2 not in O:
            O[n2] = []
        O[n2].append((t,n1))
    return O
        
def getQ(timestamps):
    O = {}
    for (t,n1,n2) in timestamps:
        if n1 not in O:
            O[n1] = []
        O[n1].append(t)
        if n2 not in O:
            O[n2] = []
        O[n2].append(t)
    O = {n: sorted(set(l)) for n, l in O.iteritems()}
    return O
    
def runBudget_fixed_budget(timestamps, b = {}):
    nodeEdgeIndex = utils.getIndex(timestamps, 'budget')
    if not b:
        b = {n: (timestamps[-1][0]-timestamps[0][0]) for n in nodeEdgeIndex}
    Xstart, Xend = budgetAlgorithm(timestamps, nodeEdgeIndex, b)
    return Xstart, Xend
    
