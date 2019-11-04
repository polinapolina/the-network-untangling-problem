from datetime import datetime, timedelta
import numpy as np
import random
import utils_k as utils
import pickle

def runKInner(timestamps, k, maxiter = 10):
    """
    Implements Baseline algorithm

    Parameters
    ----------
    timestamps : list of tuples
        Sorted list of interactions [(t1, n1, n2), (t2, n3, n4),...], where t is timestamp, n1 and n2 are interactiong nodes.
        Nodes in the interactions are sorted lexicographically.
    k : int 
        number of activity intervals
    maxiter : int
        maximum number of interactions

    Returns
    -------
    tuple of dicts
        (dict1, dict2): dict1 of a shape {n1: {0: t10, 1: t11, 2: t12}, n2: {0: t20, 1: t21, 2: t22}, ..}, tij is a starting point of jth activity interval of node ni.
                        dict2 of a shape {n1: {0: t10, 1: t11, 2: t12}, n2: {0: t20, 1: t21, 2: t22}, ..}, tij is an ending point of jth activity interval of node ni.

    """
    nodeEdgeIndex = utils.getIndex(timestamps, 'inner')
    m = getKInitial(nodeEdgeIndex, k, timestamps[0][0], timestamps[-1][0])
    Xstart, Xend = solveKInner(timestamps, nodeEdgeIndex, m)
    
    best = utils.getCost(Xstart, Xend)
    bestSol = (Xstart, Xend)    
    
    for k in xrange(1, maxiter):
        m = getKNextInput(nodeEdgeIndex, Xstart, Xend)
        Xstart, Xend = solveKInner(timestamps, nodeEdgeIndex, m)
        
        t = utils.getCost(Xstart, Xend)
        if t < best:
            best = t 
            bestSol = (Xstart, Xend)
        else:
            break
            
    return bestSol[0], bestSol[1]


def getKInitial(nodeEdgeIndex, k, st, fin):
    m = {}
    for n, t in nodeEdgeIndex.iteritems():
        
        tstmps = sorted([tst for (tst, n1, n2) in t])
        
        gaps = [t[i+1][0] - t[i][0] for i in xrange(len(t)-1)]
        indexes = sorted(np.argsort(gaps)[-(k-1):])
        ints = []
        s = st
        for i in indexes:
            ints.append((s, t[i][0]))
            s = t[i+1][0]
        ints.append((s, fin))
        m[n] = []
        for (s,f) in ints:
            mpid_ = (np.abs(np.array(tstmps) - ((s+f)/2.0))).argmin() 
            m[n].append(tstmps[mpid_])     
    return m

    
def getKInitialTrue(nodeEdgeIndex, active_truth):
    m = {}
    for n, ints in active_truth.iteritems():
        tstmps = sorted([tst for (tst, n1, n2) in nodeEdgeIndex[n]])
        m[n] = []
        for s,f in ints:            
            mpid_ = (np.abs(np.array(tstmps) - ((s+f)/2.0))).argmin() 
            m[n].append(tstmps[mpid_])
    return m

    
def getKNextInput(nodeEdgeIndex, Xstart, Xend):
    m = {}
    for n, t in nodeEdgeIndex.iteritems():
        m[n] = []
        for i in Xstart[n].keys():
            tstmps = [st for (st, n1, n2) in t if st >= Xstart[n][i] and st <= Xend[n][i]]
            mpid_ = (np.abs(np.array(tstmps) - ((min(tstmps) + max(tstmps))/2.0))).argmin() 
            m[n].append(tstmps[mpid_])
    return m
    
    
def getKAlphas(timestamps, m):
    a = {node: [0.0]*(len(v)+1) for node, v in m.iteritems()}
    b = {node: [np.inf]*(len(v)+1) for node, v in m.iteritems()}
    
    Xstart = {node: {} for node in m.keys()}
    Xend = {node: {} for node in m.keys()}

    alphaNodes, theta, eta = {}, {}, {}
    counter = {i: 0 for i in m.keys()}
    for (t, u, v) in timestamps:    
        
        if counter[u] < len(m[u]) and m[u][counter[u]] < t:
            counter[u] += 1
        if counter[v] < len(m[v]) and m[v][counter[v]] < t:
            counter[v] += 1            
        
        prevMu = m[u][counter[u]-1] if counter[u] > 0 else -np.Inf        
        nextMu = m[u][counter[u]] if counter[u] < len(m[u]) else np.Inf
        thetau = t - prevMu
        etau = nextMu - t        
        z1 = min([thetau - a[u][counter[u]], etau, b[u][counter[u]]])        
        
        prevMv = m[v][counter[v]-1] if counter[v] > 0 else -np.Inf
        nextMv = m[v][counter[v]] if counter[v] < len(m[v]) else np.Inf
        thetav = t - prevMv
        etav = nextMv - t        
        z2 = min([thetav - a[v][counter[v]], etav, b[v][counter[v]]])
        
                
        avar = min(z1, z2)
        
        if u not in alphaNodes:
            alphaNodes[u] = {}
            theta[u] = {}
            eta[u] = {}
        if t not in alphaNodes[u]:
            alphaNodes[u][t] = 0.0
            theta[u][t] = 0.0
            eta[u][t] = 0.0
        alphaNodes[u][t] += avar
        theta[u][t] = thetau
        eta[u][t] = etau
        
        if v not in alphaNodes:
            alphaNodes[v] = {}
            theta[v] = {}
            eta[v] = {}
        if t not in alphaNodes[v]:
            alphaNodes[v][t] = 0.0
            theta[v][t] = 0.0
            eta[v][t] = 0.0
        alphaNodes[v][t] += avar
        theta[v][t] = thetav
        eta[v][t] = etav
        
        
        a[u][counter[u]] = a[u][counter[u]] + avar
        a[v][counter[v]] = a[v][counter[v]] + avar
        
        b[u][counter[u]] = min(b[u][counter[u]] - avar, etau - avar)
        b[v][counter[v]] = min(b[v][counter[v]] - avar, etav - avar)
        
        if counter[u] < len(m[u]):
            if b[u][counter[u]] == 0 and counter[u] not in Xstart[u]:
                Xstart[u][counter[u]] = t   
        
        if counter[u] < len(m[u]) and t == m[u][counter[u]]:
            Xend[u][counter[u]] = t       
        if counter[u] > 0:
            if a[u][counter[u]] == thetau and (counter[u] not in Xstart[u] or counter[u]-1 not in Xend[u]):
                Xend[u][counter[u]-1] = t
        
        if counter[v] < len(m[v]):
            if b[v][counter[v]] == 0 and counter[v] not in Xstart[v]:
                Xstart[v][counter[v]] = t
        
        if counter[v] < len(m[v]) and t == m[v][counter[v]]:
            Xend[v][counter[v]] = t
        if counter[v] > 0:
            if a[v][counter[v]] == thetav and (counter[v] not in Xstart[v] or counter[v]-1 not in Xend[v]):
                Xend[v][counter[v]-1] = t

    return Xstart, Xend
    
def solveKInner(timestamps, nodeEdgeIndex, m):
    Xstart, Xend = getKAlphas(timestamps, m)
    
    return Xstart, Xend
    
    
if __name__ == "__main__":

    # dataset = 'random'
    # event_length = 100
    # overlap = 0.5
    # num_nodes = 100
    
    dataset = 'random'
    event_length = 10
    overlap = 1.0
    num_nodes = 10
    active_percent = 1.0
    k = 10
    
    
    #print 'set event length:', event_length
    #print 'set event overlap:', overlap
    
    G = utils.generateGraph(n = num_nodes)
    #G = utils.generateToyGraph()
    ##print 'number of nodes in the background network:', G.number_of_nodes()
    ##print 'number of edges in the background network:', G.number_of_edges()   
       
    timestamps, active_truth = utils.generateIntervals(G, event_length = event_length, overlap = overlap, active_percent = active_percent, number_intervals = k)
    #pickle.dump((G, timestamps, active_truth), open( "save.p", "wb" ))
    #G, timestamps, active_truth = pickle.load(open( "save.p", "rb" ))
    #timestamps = [(0, 2, 0), (1, 2, 0), (2, 0, 2), (3, 0, 2), (4, 2, 0), (5, 2, 0), (6, 0, 1), (7, 0, 2), (8, 1, 0), (9, 1, 0), (10, 1, 0), (11, 1, 0)]
    #timestamps = [(0, 2, 0), (1, 2, 0), (2, 2, 0), (3, 2, 0)]
    
    
    ##print timestamps
    #for i in timestamps:
        #print 'time', i[0], ':', i[1], i[2]
    #for n in active_truth.keys():
        #print 'node', n, ':', active_truth[n]

    ##print 'number of timestamps', len(timestamps)
    #timestamps = [(0, 3, 1), (1, 3, 1), (2, 1, 2), (3, 1, 3), (4, 3, 1), (5, 3, 1), (6, 1, 3), (7, 1, 2), (8, 2, 1), (9, 2, 1), (10, 2, 1), (11, 2, 1)]
    #active_truth = {1: [(2, 3), (6, 7)], 2: [(8, 9), (10, 11)], 3: [(0, 1), (4, 5)]}
    
    #timestamps  = [(1, 3, 1), (2, 1, 2), (3, 1, 2), (4, 3, 1), (5, 3, 1), (6, 1, 2), (7, 1, 2), (8, 2, 1)]
        
    #nodeEdgeIndex = utils.getIndex(timestamps, 'inner')
    ##print nodeEdgeIndex[1]
    #m = getKInitialTrue(nodeEdgeIndex, active_truth)
    ##print m
    #m = getKInitial(nodeEdgeIndex, k, timestamps[0][0], timestamps[-1][0]) 
    #print m

    #exit()    
    Xstart, Xend = runKInner(timestamps, 10)
    #Xstart, Xend = iterateKInner(timestamps, nodeEdgeIndex, k)
    #for n in Xstart.keys():
        #print n, ':',
        #print Xstart[n]
        #print  Xend[n]
        # for i in xrange(len(Xstart[n])):
            # #print Xstart[n][i], Xend[n][i], ';',
        # #print
    ##print 'start', Xstart
    # #print 'end',  Xend
    
    ##print 'true:', active_truth
    ##print 'm:', m
    uncov = utils.checkCoverage(Xstart, Xend, timestamps)
    #print uncov
    #print 'uncovered', len(uncov)
    
    ##print 'relative total length of solution =', utils.getCost(Xstart, Xend)/((event_length-1)*num_nodes)
    #print 'total length of solution:', utils.getCost(Xstart, Xend), ((event_length-1)*num_nodes*k)
    #print 'relative maximum length of solution =', utils.getMax(Xstart, Xend)/(event_length-1)
    p, r, f = utils.compareGT(Xstart,  Xend, active_truth, timestamps)
    #print 'precision =', p
    #print 'recall =', r
    #print 'f-measure =', f
    print p,r,f
    



