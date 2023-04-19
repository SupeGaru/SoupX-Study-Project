# The given R code has three functions:

# expandClusters: expands soup counts calculated at the cluster level to the cell level. It takes as input a matrix 
# of genes (rows) by clusters (columns) where counts are the number of soup counts for that gene/cluster combination 
# (clustSoupCnts), a matrix of genes (rows) by cells (columns) giving the observed counts (cellObsCnts), a mapping from 
# cells to clusters (clusters), and a weighting to give to each cell when distributing counts (cellWeights). The function 
# then determines a most likely allocation of soup counts at the cell level. The function returns a matrix of genes (rows) by 
# cells (columns) giving the number of soup counts estimated for each cell.

# initProgBar: creates a progress bar that won't ruin log files and shows progress towards 100%. It takes the minimum and
#  maximum values of a parameter and returns a txtProgressBar object to use updating progress.

# alloc: allocates tgt of something to length(bucketLims) different "buckets" subject to the constraint that each bucket
#  has a maximum value of bucketLims that cannot be exceeded. By default, counts are distributed equally between buckets,
#   but weights can be provided using ws to have the redistribution prefer certain buckets over others. The function returns
#    a vector of the same length as bucketLims containing values distributed into buckets.

import numpy as np
from scipy.sparse import dok_matrix, csr_matrix

def expand_clusters(clustSoupCnts, cellObsCnts, clusters, cellWeights, verbose=1):
    ws = cellWeights
    if verbose > 0:
        print(f"Expanding counts from {clustSoupCnts.shape[1]} clusters to {cellObsCnts.shape[1]} cells.")
    out = []
    for j in range(clustSoupCnts.shape[1]):
        if verbose > 1:
            print(f"Expanding cluster {j}")
        wCells = np.where(clusters == clustSoupCnts.columns[j])[0]
        ww = ws[wCells] / np.sum(ws[wCells])
        lims = cellObsCnts[:, wCells]
        nSoup = clustSoupCnts.iloc[:, j]
        expCnts = dok_matrix(lims.shape, dtype=np.float64)
        expCnts._update(np.array(lims.nonzero()).T, np.array(lims[lims.nonzero()])[0])
        wGenes = np.where((nSoup > 0) & (nSoup < np.sum(lims, axis=1)))[0]
        w = np.where(np.isin(expCnts.row, wGenes))[0]
        rows = np.split(w, np.unique(expCnts.row[w], return_index=True)[1][1:])
        tmp = []
        for row in rows:
            wCells = expCnts.col[row]
            weights = ww[wCells]
            n = nSoup[wCells[0]]
            counts = np.array(expCnts[row, wCells])[0]
            tmp.append(alloc(n, counts, weights))
        expCnts[expCnts.row.isin(wGenes), expCnts.col.isin(wCells)] = np.concatenate(tmp)
        out.append(csr_matrix(expCnts))
    out = np.column_stack(out)
    out = out[:, cellObsCnts.columns]
    return out

def init_prog_bar(min_val, max_val):
    print('0%   10   20   30   40   50   60   70   80   90   100%')
    print('|----|----|----|----|----|----|----|----|----|----|')
    return range(min_val, max_val)

def alloc(tgt, bucketLims, ws=None):
    """
    Allocate values to "buckets" subject to weights and constraints
    
    Allocates tgt of something to len(bucketLims) different "buckets" subject to the constraint that each bucket 
    has a maximum value of bucketLims that cannot be exceeded. By default counts are distributed equally between 
    buckets, but weights can be provided using ws to have the redistribution prefer certain buckets over others.
    
    Args:
    tgt (float): Value to distribute between buckets.
    bucketLims (list): The maximum value that each bucket can take. Must be a list of positive values.
    ws (list, optional): Weights to be used for each bucket. Default value makes all buckets equally likely.
    
    Returns:
    A list of the same length as bucketLims containing values distributed into buckets.
    """
    
    # Normalize weights
    if ws is None:
        ws = [1/len(bucketLims)] * len(bucketLims)
    else:
        ws = [w/sum(ws) for w in ws]
    
    # Save time in line
    if all(tgt * np.array(ws) <= np.array(bucketLims)):
        return [tgt * w for w in ws]
    
    # Need to order things in the order they'll be removed as the tgt increases
    o = sorted(range(len(bucketLims)), key=lambda i: bucketLims[i] / ws[i])
    w = [ws[i] for i in o]
    y = [bucketLims[i] for i in o]
    
    # The formula for number removed at entry i is
    # k_i = y_i / w_i * (1 - sum_j=0^(i-1) w_j) + sum_j=0^(i-1) y_j
    cw = np.cumsum([0] + w[:-1])
    cy = np.cumsum([0] + y[:-1])
    k = [y[i]/w[i] * (1 - cw[i]) + cy[i] for i in range(len(y))]
    
    # Handle zero-weights appropriately
    k = [float('inf') if w[i]==0 else k[i] for i in range(len(w))]
    
    # Everything that has k<=tgt will be set to y
    b = [k[i] <= tgt for i in range(len(k))]
    
    # We then need to work out how many counts to distribute we have left over and distribute them according to re-normalised weights
    resid = tgt - sum([y[i] for i in range(len(y)) if b[i]])
    w = [w[i]/(1 - sum([w[j] for j in range(len(w)) if b[j]])) for i in range(len(w))]
    out = [y[i] if b[i] else resid * w[i] for i in range(len(y))]
    
    # Need to reverse sort
    return [out[o.index(i)] for i in range(len(o))]
