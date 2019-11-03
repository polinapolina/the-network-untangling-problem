from datetime import datetime, timedelta
import numpy as np
import matplotlib.pyplot as plt
import utils
import random
import budget as budget 
import pickle
import numpy as np
import inner_point as inner
import baseline
import sys
      
if __name__ == "__main__":

    plt.rcParams.update({'font.size': 20, 'lines.linewidth': 3})
    plt.rcParams['xtick.labelsize'] = 20
    plt.rcParams['ytick.labelsize'] = 20

    event_length = int(sys.argv[1]) if len(sys.argv) >= 2  else 100
    num_nodes = int(sys.argv[2]) if len(sys.argv) >= 3  else 100

    dataset = 'random'
    
    overlaps = np.linspace(0.5, 1, 30)
           
    
    G = utils.generateGraph(n = num_nodes)
    
    Cost_IP, P_IP, R_IP, F_IP, Costmax_IP = [], [], [], [], []
    Cost_B, P_B, R_B, F_B, Costmax_B = [], [], [], [], []
    Cost_BL, P_BL, R_BL, F_BL, Costmax_BL = [], [], [], [], []
    
    for overlap in overlaps:    
        timestamps, active_truth = utils.generateIntervals(G, event_length = event_length, overlap = overlap)        
        
        nodeEdgeIndex = utils.getIndex(timestamps, 'inner')
        Xstart, Xend = inner.iterateInner(timestamps, nodeEdgeIndex)
        
        Cost_IP.append(utils.getCost(Xstart, Xend)/((event_length-1)*num_nodes))
        p, r, f = utils.compareGT(Xstart,  Xend, active_truth, timestamps)
        P_IP.append(p)
        R_IP.append(r)
        F_IP.append(f)
        Costmax_IP.append(utils.getMax(Xstart, Xend)/(event_length-1))
        
        nodeEdgeIndex = utils.getIndex(timestamps, 'budget')
        Xstart, Xend = budget.runBudget(timestamps, nodeEdgeIndex, b = {node: event_length-1 for node in active_truth})
        
        Cost_B.append(utils.getCost(Xstart, Xend)/((event_length-1)*num_nodes))
        p, r, f = utils.compareGT(Xstart,  Xend, active_truth, timestamps)
        P_B.append(p)
        R_B.append(r)
        F_B.append(f)
        Costmax_B.append(utils.getMax(Xstart, Xend)/(event_length-1))

        Xstart, Xend = baseline.baseline2(timestamps)

        Cost_BL.append(utils.getCost(Xstart, Xend)/((event_length-1)*num_nodes))
        p, r, f = utils.compareGT(Xstart,  Xend, active_truth, timestamps)
        P_BL.append(p)
        R_BL.append(r)
        F_BL.append(f)
        Costmax_BL.append(utils.getMax(Xstart, Xend)/(event_length-1))
        
    name = 'test_overlaps'    
    plt.figure()
    plt.plot(overlaps, F_IP, 'k')
    plt.plot(overlaps, F_B, '--r')
    plt.plot(overlaps, F_BL, ':b')
    plt.xlabel('percent of overlaps', fontsize=25)
    plt.legend(['Inner', 'Budget', '1-Baseline'], loc=0)
    plt.ylabel('F-measure', fontsize = 25)
    plt.ylim((0.3, 1.005))
    plt.xlim((min(overlaps), max(overlaps)))
    plt.tight_layout()
    plt.savefig(name+'_'+'Fmeasure.pdf')    
    
    plt.figure()    
    plt.plot(overlaps, Cost_IP, 'k')
    plt.plot(overlaps, Cost_B, '--r')
    plt.plot(overlaps, Cost_BL, ':b')
    plt.legend(['Inner', 'Budget', '1-Baseline'], loc=0)
    plt.xlabel('percent of overlaps', fontsize = 25)
    plt.ylabel('Relative total length', fontsize = 25)    
    plt.xlim((min(overlaps), max(overlaps)))
    plt.tight_layout()
    plt.savefig(name+'_'+'rel_total.pdf')
    
    plt.figure()    
    plt.plot(overlaps, Costmax_IP, 'k')
    plt.plot(overlaps, Costmax_B, '--r')
    plt.plot(overlaps, Costmax_BL, ':b')
    plt.legend(['Inner', 'Budget', '1-Baseline'], loc=0)
    plt.xlabel('percent of overlaps', fontsize=25)
    plt.ylabel('Relative max length', fontsize = 25)
    plt.xlim((min(overlaps), max(overlaps)))
    plt.ylim((0.7, 30))
    plt.yticks([1]+range(5,31,5))
    plt.tight_layout()
    plt.savefig(name+'_'+'rel_max.pdf')
    plt.show()
    

