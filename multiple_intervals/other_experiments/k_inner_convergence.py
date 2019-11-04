from datetime import datetime, timedelta
import numpy as np
import matplotlib.pyplot as plt
import copy
import random
import pickle
import argparse
import numpy as np

import sys
sys.path.append("..")
import k_inner_point as inner
import k_utils as utils


      
if __name__ == "__main__":

    plt.rcParams.update({'font.size': 20, 'lines.linewidth': 3})
    plt.rcParams['xtick.labelsize'] = 20
    plt.rcParams['ytick.labelsize'] = 20

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-k", help="number of intervals", default=10, type=int)
    parser.add_argument("--intlen", help="length of each active interval", default=100, type=int)
    parser.add_argument("--overlap", help="activity intervals overlap parameter, between 0 and 1", default=0.5, type=float)
    parser.add_argument("--nnodes", help="number of nodes in the graph", default=100, type=int)
    args = parser.parse_args()

    event_length = args.intlen
    overlap = args.overlap
    num_nodes = args.nnodes
    k = args.k

    num_of_iters = 10
     
    iterations = range(1,num_of_iters + 1)
    G = utils.generateGraph(n = num_nodes)
    
    
    timestamps, active_truth = utils.generateIntervals(G, event_length = event_length, distance_inter = 1, overlap = overlap, number_intervals = k)
    nodeEdgeIndex = utils.getIndex(timestamps, 'inner')    
        
    Cost_IP, P_IP, R_IP, F_IP, Costmax_IP = [], [], [], [], []
        
    for iter in iterations:
        if iter == 1:
            m = inner.getKInitial(nodeEdgeIndex, k, timestamps[0][0], timestamps[-1][0])        
            Xstart, Xend = inner.solveKInner(timestamps, nodeEdgeIndex, m)    
        else:
            m = inner.getKNextInput(nodeEdgeIndex, Xstart, Xend)
            Xstart, Xend = inner.solveKInner(timestamps, nodeEdgeIndex, m)       
        
        Cost_IP.append(utils.getCost(Xstart, Xend)/((event_length-1)*num_nodes*k))
        p, r, f = utils.compareGT(Xstart,  Xend, active_truth, timestamps)
        P_IP.append(p)
        R_IP.append(r)
        F_IP.append(f)
        Costmax_IP.append(utils.getMax(Xstart, Xend)/(event_length-1))
  
        
    
    name = 'test_iterations_inner'
    plt.figure()
    plt.plot(iterations, Cost_IP, 'k')
    plt.xlabel('iterations', fontsize=25)    
    plt.ylabel('Relative total length', fontsize=25)
    plt.tight_layout()

    plt.figure()
    plt.plot(iterations, Costmax_IP, 'k')
    plt.xlabel('iterations', fontsize=25)
    plt.ylabel('Relative max length', fontsize=25)
    plt.tight_layout()
    
    plt.figure()
    plt.plot(iterations, P_IP, 'b--')
    plt.plot(iterations, R_IP, 'k')
    plt.plot(iterations, F_IP, 'r:')
    plt.legend(['P', 'R', 'F'], loc = 0)
    plt.xlabel('iterations', fontsize=25)
    plt.ylabel('Quality', fontsize=25)
    plt.tight_layout()
    plt.show()
    

