from datetime import datetime, timedelta
import numpy as np
import matplotlib.pyplot as plt
import random
import pickle
import numpy as np
import copy
import uuid
import argparse

import sys
sys.path.append("..")
import k_utils as utils
import k_budget as budget 
import k_inner_point as inner


def messup_inner(m, k, rand_degree, nodeEdgeIndex):
    messed_m = {}
    all = len(m)*k
    messed_fraction = int(all*rand_degree)    
    msk = np.array([1]*messed_fraction + [0]*(all - messed_fraction))
    np.random.shuffle(msk)
    counter = 0
    for n in m.keys():
        messed_m[n] = []
        for i in xrange(len(m[n])):
            if msk[counter] == 1:
                new_val = np.random.choice([j[0] for j in nodeEdgeIndex[n]])                
            else:
                new_val = m[n][i]
            while new_val in messed_m[n]:
                new_val = np.random.choice([j[0] for j in nodeEdgeIndex[n]])
            messed_m[n].append(new_val)
            counter += 1
        messed_m[n].sort()
    return messed_m
    
def messup_gaps(m, k, rand_degree, nodeEdgeIndex):
    messed_m = {}
    all = len(m)*(k-1)
    messed_fraction = int(all*rand_degree)    
    msk = np.array([1]*messed_fraction + [0]*(all - messed_fraction))
    np.random.shuffle(msk)
    counter = 0
    for n in m.keys():
        messed_m[n] = []
        for i in xrange(len(m[n])):
            if msk[counter] == 1:
                #messed_m[n].append(np.random.randint(min_t, high = max_t+1))
                new_val = np.random.choice(nodeEdgeIndex[n].keys())
            else:
                new_val = m[n][i]
            while new_val in messed_m[n]:
                new_val = np.random.choice(nodeEdgeIndex[n].keys())
            messed_m[n].append(new_val)
            counter += 1
        messed_m[n].sort()
    return messed_m
    
def compareGT(Xstart_zero, Xend_zero, Xstart, Xend, timestamps):
    count, total = 0, 0.
    for (t, n1, n2) in timestamps:
        for n in [n1,n2]:
            total += 1
            ts, ts_zero = 0,0
            for i in xrange(len(Xstart[n])):
                if t >= Xstart[n][i] and t <= Xend[n][i]:
                    ts = 1
                    
            for i in xrange(len(Xstart[n])):
                if t >= Xstart_zero[n][i] and t <= Xend_zero[n][i]:
                    ts_zero = 1
            if ts != ts_zero:
                count += 1
    return float(count)/total 
    

      
if __name__ == "__main__":

    plt.rcParams.update({'font.size': 20, 'lines.linewidth': 3})
    plt.rcParams['xtick.labelsize'] = 20
    plt.rcParams['ytick.labelsize'] = 20

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-k", help="number of intervals", default=10, type=int)
    parser.add_argument("--intlen", help="length of each active interval", default=100, type=int)
    parser.add_argument("--overlap", help="activity intervals overlap parameter, between 0 and 1", default=0.5, type=float)
    parser.add_argument("--nnodes", help="number of nodes in the graph", default=100, type=int)
    parser.add_argument("--npoints", help="number of different values of distortion between 0 and 1 to try", default=100, type=int)
    args = parser.parse_args()

    k = args.k
    event_length = args.intlen
    overlap = args.overlap
    num_nodes = args.nnodes

    randomness = np.linspace(0.0, 1.0, args.npoints)
    
    G = utils.generateGraph(n = num_nodes)
    
    Cost_IP, F_IP, Hdist_IP = [], [], []
    Cost_B, F_B, Hdist_B = [], [], []
    
    timestamps, active_truth = utils.generateIntervals(G, event_length = event_length, overlap = overlap, number_intervals = k)
    nodeEdgeIndex = utils.getIndex(timestamps, 'inner')
    m = inner.getKInitialTrue(nodeEdgeIndex, active_truth)
    XstartI, XendI = inner.solveKInner(timestamps, nodeEdgeIndex, m)
    
    nodeEdgeIndex = utils.getIndex(timestamps, 'budget')
    m = budget.getTrueGapPoints(nodeEdgeIndex,  active_truth)
    XstartB, XendB = budget.runKBudget_BS(timestamps,nodeEdgeIndex, k, maxiters = 10, m = m)
    
    for init_rand in randomness[1:]:
        
        nodeEdgeIndex = utils.getIndex(timestamps, 'inner')
        m = inner.getKInitialTrue(nodeEdgeIndex, active_truth)
        m = messup_inner(m, k, init_rand, nodeEdgeIndex)
        
        Xstart, Xend = inner.runKInner(timestamps, k, maxiter = 10, m = m)
        
        Cost_IP.append(float(utils.getCost(Xstart, Xend))/((event_length-1)*num_nodes*k))
        p, r, f = utils.compareGT(Xstart,  Xend, active_truth, timestamps)
        F_IP.append(f)
        Hdist_IP.append(compareGT(XstartI, XendI, Xstart, Xend, timestamps))
        
        nodeEdgeIndex = utils.getIndex(timestamps, 'budget')
        m = budget.getTrueGapPoints(nodeEdgeIndex,  active_truth)
        m = messup_gaps(m, k, init_rand, nodeEdgeIndex)
        Xstart, Xend = budget.runKBudget(timestamps, k, maxiters = 10, m = m)
        
        Cost_B.append(float(utils.getCost(Xstart, Xend))/((event_length-1)*num_nodes*k))
        p, r, f = utils.compareGT(Xstart,  Xend, active_truth, timestamps)
        F_B.append(f)
        Hdist_B.append(compareGT(XstartB, XendB, Xstart, Xend, timestamps))
        
    
    plt.figure()
    plt.plot(randomness[1:], F_IP, 'k')
    plt.plot(randomness[1:], F_B, '--r')
    plt.xlabel('distortion', fontsize=25)
    plt.legend(['k-Inner', 'k-Budget'], loc=0)
    plt.ylabel('F-measure', fontsize = 25)
    plt.xlim((randomness[1], max(randomness)))
    plt.tight_layout()
    
    plt.figure()    
    plt.plot(randomness[1:], Cost_IP, 'k')
    plt.plot(randomness[1:], Cost_B, '--r')
    plt.legend(['k-Inner', 'k-Budget'], loc=0)
    plt.xlabel('distortion', fontsize = 25)
    plt.ylabel('Relative total length', fontsize = 25)    
    plt.xlim((randomness[1], max(randomness)))
    plt.tight_layout()
    
    plt.figure()    
    plt.plot(randomness[1:], Hdist_IP, 'k')
    plt.plot(randomness[1:], Hdist_B, '--r')
    plt.legend(['k-Inner', 'k-Budget'], loc=0)
    plt.xlabel('distortion', fontsize=25)
    plt.ylabel('Hamming distance', fontsize = 25)
    plt.xlim((randomness[1], max(randomness)))
    plt.tight_layout()
    plt.show()
    

