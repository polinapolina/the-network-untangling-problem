from datetime import datetime, timedelta
import numpy as np
import random
import utils

def getInitial(nodeEdgeIndex):
    m = {}
    for n, t in nodeEdgeIndex.iteritems():
        
        tstmps = [st for (st, n1, n2) in t]
        mn = 1.0*sum(tstmps)/len(tstmps)
        diff = [np.abs(i-mn) for i in tstmps]
        mindiff = min(diff)
        indixes = np.where(np.array(diff) == mindiff)[0]
        idx = indixes[len(indixes)/2]
        m[n] = tstmps[idx]
    return m
    
def getInitialTrue(nodeEdgeIndex, active_truth):
    m = {}
    for n, t in nodeEdgeIndex.iteritems():
        ts = [i for i,j,k in t if i >= active_truth[n][0] and i <= active_truth[n][1]]        
        m[n] = random.sample(ts,1)[0]
    return m
    
    
def getNextInput(nodeEdgeIndex, Xstart, Xend):
    m = {}
    for n, t in nodeEdgeIndex.iteritems():
        tstmps = [st for (st, n1, n2) in t if st >= Xstart[n] and st <= Xend[n]]
        mn = 1.0*sum(tstmps)/len(tstmps)
        diff = [np.abs(i-mn) for i in tstmps]
        mindiff = min(diff)
        indixes = np.where(np.array(diff) == mindiff)[0]
        idx = indixes[len(indixes)/2]
        m[n] = tstmps[idx]
    return m
    
    
def getAlphas(timestamps, m):
    a = {i: 0.0 for i in m}
    b = {i: np.inf for i in m}
    #tstmps = sorted(timeSlackIndex)
    alpha = {}
    alphaNodes = {}
    for (t, u, v) in timestamps:
        if t < m[v]:
            x = min(m[v] - t, b[v])
        else:
            x = t - m[v] - a[v]
            
        if t < m[u]:
            y = min(m[u] - t, b[u])
        else:
            y = t - m[u] - a[u]            
        
        delta = min(x,y)
        alpha[(t, u, v)] = delta
        
        if u not in alphaNodes:
            alphaNodes[u] = {}
        if t not in alphaNodes[u]:
            alphaNodes[u][t] = 0.0
        alphaNodes[u][t] += delta
        
        if v not in alphaNodes:
            alphaNodes[v] = {}
        if t not in alphaNodes[v]:
            alphaNodes[v][t] = 0.0
        alphaNodes[v][t] += delta
        
        if t <  m[v]:   
            b[v] = min(b[v] - delta, m[v] - t - delta)
        else:
            a[v] = a[v] + delta
            
        if t <  m[u]:
            b[u] = min(b[u] - delta, m[u] - t - delta)
        else:
            a[u] = a[u] + delta
        
    return alpha, alphaNodes
    
def getSolution(alphaNodes, m):
    Xstart = {}
    Xend = {}
    
    for n, ts in alphaNodes.iteritems():
        k = sorted(ts.keys())
        vals = [ts[i] for i in k]
        ind = k.index(m[n])
        for i in xrange(0, ind+1):
            if sum(vals[i:ind+1]) == (m[n] - k[i]):
                Xstart[n] = k[i]
                break
                
        for i in xrange(len(k)-1, ind-1, -1):
            if sum(vals[ind:i+1]) == (k[i] - m[n]):
                Xend[n] = k[i]
                break
        if n not in Xend:
            Xend[n] = m[n]
        if n not in Xstart:
            Xstart[n] = m[n]
    return Xstart, Xend
 
    
def solveInner(timestamps, nodeEdgeIndex, m):
    alpha, alphaNodes = getAlphas(timestamps, m)
    Xstart, Xend = getSolution(alphaNodes, m)
    return Xstart, Xend
    
def iterateInner(timestamps, nodeEdgeIndex, max_iter = 10):

    m = getInitial(nodeEdgeIndex)
    alpha, alphaNodes = getAlphas(timestamps, m)
    Xstart, Xend = getSolution(alphaNodes, m)
    
    best = utils.getCost(Xstart, Xend)
    bestSol = (Xstart, Xend)    
    
    for k in xrange(1, max_iter):
        m = getNextInput(nodeEdgeIndex, Xstart, Xend)
        alpha, alphaNodes = getAlphas(timestamps, m)
        Xstart, Xend = getSolution(alphaNodes, m)
        
        t = utils.getCost(Xstart, Xend)
        if t < best:
            best = t 
            bestSol = (Xstart, Xend)
        else:
            break
            
    return bestSol
    
    
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
        
    nodeEdgeIndex = utils.getIndex(timestamps, 'inner')    
    m = getInitialTrue(nodeEdgeIndex, active_truth)    
    Xstart, Xend = solveInner(timestamps, nodeEdgeIndex, m)
    print
    
    print 'relative total length of solution =', utils.getCost(Xstart, Xend)/((event_length-1)*num_nodes)
    print 'relative maximum length of solution =', utils.getMax(Xstart, Xend)/(event_length-1)
    p, r, f = utils.compareGT(Xstart,  Xend, active_truth, timestamps)
    print 'precision =', p
    print 'recall =', r
    print 'f-measure =', f
    



