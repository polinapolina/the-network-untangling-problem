from datetime import datetime, timedelta
import numpy as np
import matplotlib.pyplot as plt
import utils
import random
import pickle
import numpy as np
import sys
import inner_point as inner


      
if __name__ == "__main__":

    plt.rcParams.update({'font.size': 20, 'lines.linewidth': 3})
    plt.rcParams['xtick.labelsize'] = 20
    plt.rcParams['ytick.labelsize'] = 20

    dataset = 'random'

    event_length = int(sys.argv[1]) if len(sys.argv) >= 2 else 100
    overlap = float(sys.argv[2]) if len(sys.argv) >= 3 else 0.5
    num_nodes = int(sys.argv[3]) if len(sys.argv) >= 4  else 100
    
    
    G = utils.generateGraph(n = num_nodes)
    timestamps, active_truth = utils.generateIntervals(G, event_length = event_length, overlap = overlap)
    nodeEdgeIndex = utils.getIndex(timestamps, 'inner')    
    
    Cost_IP, P_IP, R_IP, F_IP, Costmax_IP = [], [], [], [], []
    
    iterations = range(1,11)
    
    for iter in iterations:
        if iter == 1:
            m = inner.getInitial(nodeEdgeIndex)
            Xstart, Xend = inner.solveInner(timestamps, nodeEdgeIndex, m)
        else:
            m = inner.getNextInput(nodeEdgeIndex, Xstart, Xend)
            Xstart, Xend = inner.solveInner(timestamps, nodeEdgeIndex, m)        
        
        Cost_IP.append(utils.getCost(Xstart, Xend)/((event_length-1)*num_nodes))
        p, r, f = utils.compareGT(Xstart,  Xend, active_truth, timestamps)
        P_IP.append(p)
        R_IP.append(r)
        F_IP.append(f)
        Costmax_IP.append(utils.getMax(Xstart, Xend)/(event_length-1))
        
    name = 'test_convergence'
    plt.figure()
    plt.plot(iterations, Cost_IP, 'k')
    plt.xlabel('iterations', fontsize=25)    
    plt.ylabel('Relative total length', fontsize=25)
    plt.tight_layout()
    plt.savefig(name+'_'+'total.pdf')

    plt.figure()
    plt.plot(iterations,Costmax_IP, 'k')
    plt.xlabel('iterations', fontsize=25)
    plt.ylabel('Relative max length', fontsize=25)
    plt.tight_layout()
    plt.savefig(name+'_'+'costmax.pdf')
    
    plt.figure()
    plt.plot(iterations, P_IP, 'b--')
    plt.plot(iterations, R_IP, 'k')
    plt.plot(iterations, F_IP, 'r:')
    plt.legend(['P', 'R', 'F'], loc = 0)
    plt.xlabel('iterations', fontsize=25)
    plt.ylabel('Quality', fontsize=25)
    plt.tight_layout()
    plt.savefig(name+'_'+'quality.pdf')
    plt.show()
    

