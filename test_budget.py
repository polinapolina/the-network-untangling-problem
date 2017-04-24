from datetime import datetime, timedelta
import numpy as np
import copy
import sys
import pickle
import utils

sys.setrecursionlimit(10000000)


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
                DFS(n, t, nodeEdgeIndex, O,  b, Smain)
    clearVisited(nodeEdgeIndex)
    Q = getQ(timestamps)
    R = []
    for (node, tm) in Smain[::-1]:
        if not nodeEdgeIndex[node][tm]['visited']:
            reverseDFS(node, tm, nodeEdgeIndex, Q, b ,R)
                
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
    
def runBudget(timestamps, nodeEdgeIndex, b = {}):
    #nodeEdgeIndex = utils.getIndex(timestamps, alg='budget')
    if not b:
        b = {n: (timestamps[-1][0]-timestamps[0][0]) for n in nodeEdgeIndex}
    Xstart, Xend = budgetAlgorithm(timestamps, nodeEdgeIndex, b)
    return Xstart, Xend
    
def runBudget_BS(timestamps,nodeEdgeIndex, klow = -1, kup = -1):
    
    #nodeEdgeIndex = getIndex(timestamps)
    maxb = timestamps[-1][0]-timestamps[0][0]
    if klow ==  -1:
        klow, kup = 1./maxb, 1.0 
    b = {n: maxb for n in nodeEdgeIndex}
    LXstart, LXend = budgetAlgorithm(timestamps, nodeEdgeIndex, b)
    
    alpha = 1.0/maxb
    while (kup - klow) > alpha:
        k = (kup + klow)/2
        b = {n: maxb*k for n in nodeEdgeIndex}
        Xstart, Xend = budgetAlgorithm(timestamps, nodeEdgeIndex, b)
        uncovered =  utils.checkCoverage(Xstart, Xend, timestamps)
        if uncovered:
            klow = k
        else:
            kup = k
            LXstart, LXend = Xstart, Xend
            
    return LXstart, LXend

if __name__ == "__main__":

    dataset = 'random'
    event_length = 100
    overlap = 0.5
    num_nodes = 100
    
    print 'set event length:', event_length
    print 'set event overlap:', overlap
    
    G = utils.generateGraph(n = num_nodes)
    print 'number of nodes in the background network:', G.number_of_nodes()
    print 'number of edges in the background network:', G.number_of_edges()   
       
    timestamps, active_truth = utils.generateIntervals(G, event_length = event_length, overlap = overlap)

    print 'number of timestamps', len(timestamps)
        
    nodeEdgeIndex = utils.getIndex(timestamps, 'budget')
    Xstart, Xend = runBudget(timestamps, nodeEdgeIndex, b = {node: event_length-1 for node in active_truth})
    print
    
    print 'relative total length of solution =', utils.getCost(Xstart, Xend)/((event_length-1)*num_nodes)
    print 'relative maximum length of solution =', utils.getMax(Xstart, Xend)/(event_length-1)
    p, r, f = utils.compareGT(Xstart,  Xend, active_truth, timestamps)
    print 'precision =', p
    print 'recall =', r
    print 'f-measure =', f


