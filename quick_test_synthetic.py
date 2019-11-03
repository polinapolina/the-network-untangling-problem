from datetime import datetime, timedelta
import numpy as np
import random
import sys
import pickle
import networkx as nx
import time 
import utils
import operator
import baseline, budget, inner_point
    
    
      
if __name__ == "__main__":
    dataset = 'random'

    alg = sys.argv[1]
    event_length = int(sys.argv[2]) if len(sys.argv) >= 3  else 100
    overlap = float(sys.argv[3]) if len(sys.argv) >= 4  else 0.2
    num_nodes = int(sys.argv[4]) if len(sys.argv) >= 5  else 100
    
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
        nodeEdgeIndex = utils.getIndex(timestamps, 'inner')    
        m = inner_point.getInitialTrue(nodeEdgeIndex, active_truth)    
        Xstart, Xend = inner_point.solveInner(timestamps, nodeEdgeIndex, m)
    elif alg == 'budget': 
        nodeEdgeIndex = utils.getIndex(timestamps, 'budget')
        Xstart, Xend = budget.runBudget(timestamps, nodeEdgeIndex, b = {node: event_length-1 for node in active_truth})
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


