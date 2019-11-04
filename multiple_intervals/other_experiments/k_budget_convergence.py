from datetime import datetime, timedelta
import numpy as np
import matplotlib.pyplot as plt
import pickle
import numpy as np
import copy
import uuid
import argparse
import random

import sys
sys.path.append("..")

import k_utils as utils
import k_budget as budget 

      
if __name__ == "__main__":

    plt.rcParams.update({'font.size': 20, 'lines.linewidth': 3})
    plt.rcParams['xtick.labelsize'] = 20
    plt.rcParams['ytick.labelsize'] = 20

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-k", help="number of intervals", default=10, type=int)
    parser.add_argument("--intlen", help="length of each active interval", default=10, type=int)
    parser.add_argument("--overlap", help="activity intervals overlap parameter, between 0 and 1", default=0.5, type=float)
    parser.add_argument("--nnodes", help="number of nodes in the graph", default=10, type=int)
    args = parser.parse_args()

    event_length = args.intlen
    overlap = args.overlap
    num_nodes = args.nnodes
    k = args.k
           
    G = utils.generateGraph(n = num_nodes)
    
    timestamps, active_truth = utils.generateIntervals(G, event_length = event_length, distance_inter = 1, overlap = overlap, number_intervals = k)   
    nodeEdgeIndex = utils.getIndex(timestamps, 'budget')
    
    m = budget.getGapPoints(nodeEdgeIndex, timestamps[0][0],timestamps[-1][0], k)
    Xstart,  Xend = budget.runKBudget_BS(timestamps, nodeEdgeIndex, k, maxiters = 100, m = m)
    
    Cost_B, P_B, R_B, F_B, Costmax_B = [], [], [], [], []
    Cost_B.append(utils.getCost(Xstart, Xend)/((event_length-1)*num_nodes*k))
    p, r, f = utils.compareGT(Xstart,  Xend, active_truth, timestamps)
    P_B.append(p)
    R_B.append(r)
    F_B.append(f)
    Costmax_B.append(utils.getMax(Xstart, Xend)/(event_length-1))
    iterations = range(1,11)
    
    for iter in iterations[1::]:
        m = budget.getNewGapPoints(nodeEdgeIndex, timestamps[0][0],timestamps[-1][0], Xstart,  Xend)
        Xstart,  Xend = budget.runKBudget_BS(timestamps, nodeEdgeIndex, k, maxiters = 100, m = m)
     
        Cost_B.append(utils.getCost(Xstart, Xend)/((event_length-1)*num_nodes*k))
        p, r, f = utils.compareGT(Xstart,  Xend, active_truth, timestamps)
        P_B.append(p)
        R_B.append(r)
        F_B.append(f)
        Costmax_B.append(utils.getMax(Xstart, Xend)/(event_length-1))
    
    plt.figure()    
    plt.plot(iterations, Cost_B, 'k')
    plt.xlabel('iteration', fontsize = 25)
    plt.ylabel('Relative total length', fontsize = 25)
    plt.tight_layout()
    
    plt.figure()    
    plt.plot(iterations, Costmax_B, 'k')
    plt.xlabel('iteration', fontsize=25)
    plt.ylabel('Relative max length', fontsize = 25)
    plt.tight_layout()
    
    plt.figure()
    plt.plot(iterations, P_B, 'b--')
    plt.plot(iterations, R_B, 'k')
    plt.plot(iterations, F_B, 'r:')
    plt.legend(['P', 'R', 'F'], loc = 0)
    plt.xlabel('iterations', fontsize=25)
    plt.ylabel('Quality', fontsize=25)
    plt.tight_layout()
    plt.show()

