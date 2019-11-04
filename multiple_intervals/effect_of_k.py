from datetime import datetime, timedelta
import numpy as np
import matplotlib.pyplot as plt
import k_utils as utils
import random
import pickle
import numpy as np
import copy
import uuid
import argparse
import k_baseline, k_budget, k_inner_point

      
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

    k = args.k
    event_length = args.intlen
    overlap = args.overlap
    num_nodes = args.nnodes

  
    krange = range(2,11)
    
    G = utils.generateGraph(n = num_nodes)
    
    Cost_IP, P_IP, R_IP, F_IP, Costmax_IP = [], [], [], [], []
    Cost_B, P_B, R_B, F_B, Costmax_B = [], [], [], [], []
    Cost_bl, P_bl, R_bl, F_bl, Costmax_bl = [], [], [], [], []
    
    
    for k in krange:
    
        timestamps, active_truth = utils.generateIntervals(G, event_length = event_length, overlap = overlap, number_intervals = k)
        
        Xstart, Xend = k_inner_point.runKInner(timestamps,k)
        
        Cost_IP.append(float(utils.getCost(Xstart, Xend))/((event_length-1)*num_nodes*k))
        p, r, f = utils.compareGT(Xstart,  Xend, active_truth, timestamps)
        P_IP.append(p)
        R_IP.append(r)
        F_IP.append(f)
        Costmax_IP.append(float(utils.getMax(Xstart, Xend))/(event_length-1))
        
        Xstart, Xend = k_budget.runKBudget(timestamps, k)
        
        
        Cost_B.append(float(utils.getCost(Xstart, Xend))/((event_length-1)*num_nodes*k))
        p, r, f = utils.compareGT(Xstart,  Xend, active_truth, timestamps)
        P_B.append(p)
        R_B.append(r)
        F_B.append(f)
        Costmax_B.append(float(utils.getMax(Xstart, Xend))/(event_length-1))

        Xstart, Xend = k_baseline.kbaseline(timestamps, k)


        Cost_bl.append(float(utils.getCost(Xstart, Xend))/((event_length-1)*num_nodes*k))
        p, r, f = utils.compareGT(Xstart,  Xend, active_truth, timestamps)
        P_bl.append(p)
        R_bl.append(r)
        F_bl.append(f)
        Costmax_bl.append(float(utils.getMax(Xstart, Xend))/(event_length-1))
   
    
    name = 'test_krange'    
    plt.figure()
   
    plt.plot(krange, F_IP, 'k')
    plt.plot(krange, F_B, '--r')
    plt.plot(krange, F_bl, ':b')
    plt.xlabel('number of intervals', fontsize=25)
    plt.legend(['k-Inner', 'k-Budget', 'k-Baseline'], loc=0)
    plt.xlim((min(krange), max(krange)))
    plt.tight_layout()
    plt.savefig(name+'_'+'Fmeasure.pdf')    
    
    plt.figure()    
    plt.plot(krange, Cost_IP, 'k')
    plt.plot(krange, Cost_B, '--r')
    plt.plot(krange, Cost_bl, ':b')
    plt.legend(['Inner', 'Budget'], loc=0)
    plt.legend(['k-Inner', 'k-Budget', 'k-Baseline'], loc=0)
    plt.xlabel('number of intervals', fontsize=25)
    plt.ylabel('Relative total length', fontsize = 25)    
    plt.xlim((min(krange), max(krange)))
    plt.tight_layout()
    plt.savefig(name+'_'+'rel_total.pdf')
    
    plt.figure()    
    plt.plot(krange, Cost_IP, 'k')
    plt.plot(krange, Costmax_B, '--r')
    plt.plot(krange, Costmax_bl, ':b')
    plt.legend(['k-Inner', 'k-Budget', 'k-Baseline'], loc=0)
    plt.xlabel('number of intervals', fontsize=25)
    plt.ylabel('Relative max length', fontsize = 25)
    plt.xlim((min(krange), max(krange)))
    plt.tight_layout()
    plt.savefig(name+'_'+'rel_max.pdf')
    plt.show()
    

