from datetime import datetime, timedelta
import numpy as np
import random
import sys
import pickle
import networkx as nx
import time 
import utils
import operator
import argparse
import baseline, budget, inner_point
    
    
      
if __name__ == "__main__":

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("algorithm", help="baseline, budget, or inner")
    parser.add_argument("--intlen", help="length of each active interval", default=100, type=int)
    parser.add_argument("--overlap", help="activity intervals overlap parameter, between 0 and 1", default=0.2, type=float)
    parser.add_argument("--nnodes", help="number of nodes in the graph", default=100, type=int)
    args = parser.parse_args()

    alg = args.algorithm
    event_length = args.intlen
    overlap = args.overlap
    num_nodes = args.nnodes
    
    print 'algorithm:', alg
    print 'set event length:', event_length
    print 'set event overlap:', overlap
    
    G = utils.generateGraph(n = num_nodes)
    print 'number of nodes in the background network:', G.number_of_nodes()
    print 'number of edges in the background network:', G.number_of_edges()
       
    timestamps, active_truth = utils.generateIntervals(G, event_length = event_length, overlap = overlap)

    print 'number of timestamps', len(timestamps)

    if alg == 'baseline':
        Xstart, Xend = baseline.baseline2(timestamps)
    elif alg == 'inner':
        Xstart, Xend = inner_point.runInner(timestamps)
    elif alg == 'budget': 
        Xstart, Xend = budget.runBudget(timestamps)
    else: 
        print('no algorithm specified')
        exit()

    print
        
    print 'relative total length of solution =', utils.getCost(Xstart, Xend)/((event_length-1)*num_nodes)
    print 'relative maximum length of solution =', utils.getMax(Xstart, Xend)/(event_length-1)
    p, r, f = utils.compareGT(Xstart,  Xend, active_truth, timestamps)
    print 'precision =', p
    print 'recall =', r
    print 'f-measure =', f


