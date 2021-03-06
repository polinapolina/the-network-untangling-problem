from datetime import datetime, timedelta
import numpy as np
import networkx as nx
from networkx.utils import powerlaw_sequence
import random
from random import shuffle
import time

def readdata(filename, unix = True):
    timestamps = []
    with open(filename, 'r') as f:    
        for line in f:                
            line = line.strip().split(' ')
            if not unix:
                tstr =  line[0][1:] + ' ' +line[1][0:-1]
                t = datetime.strptime(tstr, '%Y-%m-%d %H:%M:%S')
                timestamp = time.mktime(datetime.strptime(tstr, '%Y-%m-%d %H:%M:%S').timetuple())
                
                tst, n1, n2 = int(timestamp), line[2], line[3]
            else:
                timestamp = int(line[0])
                tst, n1, n2 = int(line[0]), int(line[1]), int(line[2])
            
            if n1 == n2:
                pass
            
            if n2 < n1:
                n1, n2 = n2, n1
            
            timestamps.append((tst, n1, n2))
            
        timestamps.sort()
    return timestamps

    
def compareGT(Xstart, Xend, active_truth, timestamps):
    p, r = {}, {}
    TP = {}
    P = {}
    S = {}
    for (t, n1, n2) in timestamps:
        for n in [n1,n2]:
            active_ints = active_truth[n]
            #print 'n', n
            if n not in P:
                P[n] = 0.0
            if n not in S:
                S[n] = 0.0
            if n not in TP:
                TP[n] = 0.0
            
            tp = 0.0
            for (s,f) in active_ints:
                if t >= s and t <= f:
                    tp = 1.0
            P[n] += tp  
                 
            ts = 0.0
            for i in xrange(len(Xstart[n])):
                if t >= Xstart[n][i] and t <= Xend[n][i]:
                    ts = 1.0
            S[n] += ts
                
            if ts == 1.0 and tp == 1.0:
                TP[n] += 1.0             
    for n in P:
        p[n] = TP[n]/S[n] if S[n] > 0 else 0.0 
        r[n] = TP[n]/P[n] if P[n] > 0 else 0.0
    p,r = np.mean(p.values()), np.mean(r.values())
    f = (2.0*(p*r)/(p+r) if (p+r) > 0 else 0.0 )
    return p,r,f

def getIndex(timestamps, alg = 'inner'):
    nodeEdgeIndex = {}
    
    for (tst,n1, n2) in timestamps:
        if n2 < n1:
            n1, n2 = n2, n1
        
        if alg == 'inner':
            if n1 not in nodeEdgeIndex:
                nodeEdgeIndex[n1] = []
            if n2 not in nodeEdgeIndex:
                nodeEdgeIndex[n2] = []
            nodeEdgeIndex[n1].append((tst, n1, n2))
            nodeEdgeIndex[n2].append((tst, n1, n2))
        elif alg == 'budget':
            if n1 not in nodeEdgeIndex:
                nodeEdgeIndex[n1] = {}
            if n2 not in nodeEdgeIndex:
                nodeEdgeIndex[n2] = {}
                
            if tst not in nodeEdgeIndex[n1]:
                nodeEdgeIndex[n1][tst] = {'visited': False , 'edges': set() }
            nodeEdgeIndex[n1][tst]['edges'].add((n1, n2))
            if tst not in nodeEdgeIndex[n2]:
                nodeEdgeIndex[n2][tst] = {'visited': False, 'edges': set() }
            nodeEdgeIndex[n2][tst]['edges'].add((n1, n2))
            
    return nodeEdgeIndex
        

def readgraph(filename):
    timeSlackIndex, nodeEdgeIndex = {}, {}
    G = nx.Graph()
    with open(filename, 'r') as f:    
        for line in f:                
            line = line.strip().split(' ')
            if int(line[2]) != int(line[3]):
                G.add_edge(int(line[2]), int(line[3]))
            
    return G

def generateGraph(n = 100, seed = 1.0):
    generated = False
    z = nx.utils.create_degree_sequence(n, powerlaw_sequence)
    G = nx.configuration_model(z, seed = seed)
    G = nx.Graph(G)
    G.remove_edges_from(G.selfloop_edges())
    return G
    
def generateToyGraph():
    #edges = [(1,2),(2,3),(1,3)]
    edges = [(1,2),(1,3)]
    G = nx.Graph()
    G.add_edges_from(edges)
    return G
    
def generateIntervals(G, distance_in = 1, event_length = 10, distance_inter = 1, overlap = 0.1, seed = 1.0, number_intervals = 10):
    G.remove_nodes_from(nx.isolates(G))
    random.seed(seed)
    timestmps = []
    
    nodes = list(G.nodes())
    active = {n: [] for n in nodes}
    allnodes = []
    for i in range(number_intervals):
        shuffle(nodes)
        allnodes += nodes
    
    current = 0
    for n in allnodes:
        l = event_length
        s, f = current, current + distance_in*(l-1)
        randt = list(s + np.random.rand(l-2)*(f-s)) + [s, f]
        for i in xrange(l):
            neigh = list(G.neighbors(n))
            n2 = np.random.choice(neigh)
            n1, n2 = (n, n2) if n < n2 else (n2, n) 
            timestmps.append((randt[i], n1, n2))

        if overlap > 0:
            current = max(f - int((f-s)*overlap), 0)
        else:
            current = f + distance_inter
        active[n].append((s,f))
    timestmps = list(set(timestmps))
    timestmps.sort()
    return timestmps, active
    
def checkCoverage(Xstart, Xend, timeSlackIndex):
    uncovered = []    
    for (t,u,v) in timeSlackIndex:
        cov = False
        for i in xrange(len(Xstart[u])):
            if t >= Xstart[u][i] and t <= Xend[u][i]:
                cov = True
        for i in xrange(len(Xstart[v])):
            if t >= Xstart[v][i] and t <= Xend[v][i]:
                cov = True
        
        if not cov:
            uncovered.append((t,u,v))
    return uncovered
    
def getCost(Xstart, Xend):
    cost = 0.
    for n in Xstart.keys():
        for i in Xstart[n].keys():
            t = Xend[n][i] - Xstart[n][i]
            if t < np.inf and t > -np.inf:
                cost += t
    return cost
   
def getMax(Xstart, Xend):
    cost = 0.
    for n in Xstart.keys():
        for i in Xstart[n].keys():
            t = Xend[n][i] - Xstart[n][i]
            if t < np.inf and t > -np.inf:
                if t > cost:
                    cost = t
    return cost

def getMaxCover(nodeEdgeIndex):
    X = {}    
    for n, ts in nodeEdgeIndex.iteritems():
        X[n] = (max(ts.keys()) - min(ts.keys()))
    return X 
    
def getAvgRate(nodeEdgeIndex):
    X = {}    
    for n, ts in nodeEdgeIndex.iteritems():
        X[n] = 1.0*(max(ts.keys()) - min(ts.keys()))/len(ts)
    return np.mean(X.values())