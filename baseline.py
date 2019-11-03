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
    """
    Implements Budget algorithm

    Parameters
    ----------
    timestamps : list of tuples
        Sorted list of interactions [(t1, n1, n2), (t2, n3, n4),...], where t is timestamp, n1 and n2 are interactiong nodes.
        Nodes in the interactions are sorted lexicographically.
    maxiter : int
        maximum number of interactions in binary search

    Returns
    -------
    tuple of dicts
        (dict1, dict2): dict1 of a shape {n1: t1, n2: t2, ..} contains starting point t1 of activity interval of node n1, only for nodes chosen to be active.
                        dict2 of a shape {n1: t1, n2: t2, ..} contains ending point t1 of activity interval of node n1, only for nodes chosen to be active.

    """
    Xstart, Xend = {}, {}
    while timestamps:
        nodeEdgeIndex = utils.getIndex(timestamps)

        gain = {node: np.float64(len(tstmps))/(tstmps[-1][0] - tstmps[0][0]) for node, tstmps in nodeEdgeIndex.iteritems()}
        node, gain = (sorted(gain.items(), key=operator.itemgetter(1)))[-1]
        s = sorted([i[0] for i in nodeEdgeIndex[node]])
        Xstart[node], Xend[node] = s[0], s[-1]
        timestamps = [x for x in timestamps if x[1] != node and x[2] != node]

    return Xstart, Xend
