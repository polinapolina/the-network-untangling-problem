from datetime import datetime, timedelta
import numpy as np
import random
import sys
import pickle
import networkx as nx
import time 
import utils
import operator
import baseline
    
    
      
if __name__ == "__main__":
    dataset = 'random'
    event_length = 10
    overlap = 0.2
    num_nodes = 100
    
    print 'set event length:', event_length
    print 'set event overlap:', overlap
    
    G = utils.generateGraph(n = num_nodes)
    print 'number of nodes in the background network:', G.number_of_nodes()
    print 'number of edges in the background network:', G.number_of_edges()   
       
    timestamps, active_truth = utils.generateIntervals(G, event_length = event_length, overlap = overlap)

    print 'number of timestamps', len(timestamps)
    Xstart, Xend = baseline.baseline2(timestamps)
        
    print 'relative total length of solution =', utils.getCost(Xstart, Xend)/((event_length-1)*num_nodes)
    print 'relative maximum length of solution =', utils.getMax(Xstart, Xend)/(event_length-1)
    p, r, f = utils.compareGT(Xstart,  Xend, active_truth, timestamps)
    print 'precision =', p
    print 'recall =', r
    print 'f-measure =', f


