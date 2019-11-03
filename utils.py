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
            s, f = active_truth[n]
            if n not in P:
                P[n] = 0.0
            if n not in S:
                S[n] = 0.0
            if n not in TP:
                TP[n] = 0.0
            
              
            if t >= s and t <= f:
                P[n] += 1.0        
            if n in Xstart and t >= Xstart[n] and t <= Xend[n]:
                S[n] += 1.0
                
            if n in Xstart and t >= s and t <= f and t >= Xstart[n] and t <= Xend[n]:
                TP[n] += 1.0             
           
    for n in active_truth:
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
    while not generated:
        try:
            z = nx.utils.create_degree_sequence(n, powerlaw_sequence)
            G = nx.configuration_model(z, seed = seed)
            generated = True
        except:
            pass
    G = nx.Graph(G)
    G.remove_edges_from(G.selfloop_edges())
    # connectors = []
    # for nodes in nx.connected_components(G):
        # connectors.append(nodes[0])
    # if connectors:
        # G.add_path(connectors)        
    return G


    
def generateIntervals(G, distance_in = 1, event_length = 10, distance_inter = 1, overlap = 0.1, seed = 1.0):
    for i in xrange(30):
        try:
            random.seed(seed)
            timestmps = []
            active = {}
            nodes = G.nodes()
            shuffle(nodes)
             
            current = 0
            for n in nodes:        
                s = current        
                l = event_length
                for i in xrange(l):
                    neigh = G.neighbors(n)
                    timestmps.append((current, n, np.random.choice(neigh)))
                    current += distance_in
                current -= distance_in
                f = current
                if overlap > 0:
                    current = max(f - int((f-s)*overlap), 0)
                else:
                    current = f + distance_inter
                active[n] = (s,f)         
                
            timestmps.sort()
                
            return timestmps, active
        except:
            pass
    print('bad parameters for synthetic dataset')
    exit()
    
def checkCoverage(Xstart, Xend, timeSlackIndex):
    uncovered = []    
    for (t,u,v) in timeSlackIndex:
        if (t < Xstart[u] or t > Xend[u]) and (t < Xstart[v] or t > Xend[v]):
            uncovered.append((t,u,v))
    return uncovered
    
def getCost(Xstart, Xend):
    cost = sum([Xend[i] - Xstart[i] for i in Xstart])    
    return cost
   
def getMax(Xstart, Xend):
    cost = max([Xend[i] - Xstart[i] for i in Xstart])    
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