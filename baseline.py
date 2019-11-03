from datetime import datetime, timedelta
import numpy as np
import random
import sys
import pickle
import networkx as nx
import time 
import utils
import operator


def baseline2(timestamps):
    Xstart, Xend = {}, {}
    while timestamps:
        nodeEdgeIndex = utils.getIndex(timestamps)

        gain = {node: np.float64(len(tstmps))/(tstmps[-1][0] - tstmps[0][0]) for node, tstmps in nodeEdgeIndex.iteritems()}
        node, gain = (sorted(gain.items(), key=operator.itemgetter(1)))[-1]
        s = sorted([i[0] for i in nodeEdgeIndex[node]])
        Xstart[node], Xend[node] = s[0], s[-1]
        timestamps = [x for x in timestamps if x[1] != node and x[2] != node]
    return Xstart, Xend
