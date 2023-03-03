# The quickMarkers R function uses Term Frequency - Inverse Document Frequency (tf-idf) ordering to get the top N markers 
# of each cluster. For each cluster, it returns either the top N or all genes passing the hypergeometric test with the
#  False Discovery Rate (FDR) specified, whichever list is smallest.

# The function first binarizes gene expression data in each cell, replacing the counts with 1 if the count is above a
#  threshold (expressCut), and 0 otherwise. The frequency with which a gene is expressed within the target group is compared 
#  to the global frequency to calculate the tf-idf score. The function then calculates a multiple hypothesis-corrected p-value
#   based on a hypergeometric test.

# The input to the function is a table of counts (toc), a vector of cluster membership (clusters), the number of marker genes
#  to return per cluster (N), the FDR to use, and a threshold value for binarization (expressCut). The output is a data.frame
#   with top N markers (or all that pass the hypergeometric test) and their statistics for each cluster. The columns of the output
#    data.frame include gene name, cluster, gene frequency, gene frequency outside the cluster, gene frequency second best,
#    gene frequency global, second-best cluster name, tf-idf score, idf score, and q-value.

import pandas as pd
import numpy as np
from scipy.sparse import issparse
from scipy.stats import hypergeom

def quickMarkers(toc, clusters, N=10, FDR=0.01, expressCut=0.9):
    """Gets top N markers for each cluster using tf-idf ordering and hypergeometric test.

    Args:
        toc (scipy.sparse.spmatrix): Table of counts. Must be a sparse matrix.
        clusters (array-like): Vector of length n_cells giving cluster membership.
        N (int): Number of marker genes to return per cluster.
        FDR (float): False discovery rate to use.
        expressCut (float): Value above which a gene is considered expressed.

    Returns:
        pandas.DataFrame: A data frame with top N markers (or all that pass the hypergeometric test) 
                           and their statistics for each cluster.
    """
    if not issparse(toc):
        raise ValueError("toc must be a sparse matrix")

    toc = toc.tocsr()
    w = toc.data > expressCut
    gene_indices = toc.indices[w]
    cell_indices = toc.indptr[:-1][w]

    # Get the counts in each cluster
    clCnts = pd.Series(clusters).value_counts()

    gene_names = toc.indices.astype(str)
    nObs = pd.crosstab(gene_names[gene_indices], clusters[cell_indices], rownames=["Gene"], colnames=["Cluster"])
    nTot = nObs.sum(axis=1)
    tf = nObs.div(clCnts).T
    ntf = (nTot[:, np.newaxis] - nObs).div(toc.shape[1] - clCnts).T
    idf = np.log(toc.shape[1] / nTot)

    score = tf * idf
    qvals = pd.DataFrame(np.zeros_like(score), columns=score.columns, index=score.index)
    for cluster in qvals.columns:
        _, qvals[cluster], _, _ = multipletests(
            hypergeom.sf(nObs[cluster] - 1, toc.shape[1], nTot, clCnts[cluster]),
            alpha=FDR,
            method='fdr_bh'
        )

    sndBest = tf.drop(columns=score.idxmax(axis=1)).max(axis=1)
    sndBestName = tf.drop(columns=score.idxmax(axis=1)).idxmax(axis=1)

    # Now get the top N for each group
    w = []
    for col in score.columns:
        o = np.argsort(score[col])[::-1]
        if np.sum(qvals[col] < FDR) >= N:
            w.append(o[:N])
        else:
            w.append(o[qvals[col][o] < FDR])
    ww = np.c_[np.ravel(w, order='F'), np.repeat(np.arange(ncol(nObs)), np.ravel(w, order='F'))]
    out = pd.DataFrame({
    'gene': nObs.index.values[ww[:,0]],
    'cluster': nObs.columns.values[ww[:,1]],
    'geneFrequency': tf[ww[:,0], ww[:,1]],
    'geneFrequencyOutsideCluster': ntf[ww[:,0], ww[:,1]],
    'geneFrequencySecondBest': sndBest[ww[:,0], ww[:,1]],
    'geneFrequencyGlobal': nTot[ww[:,0]]/ncol(toc),
    'secondBestClusterName': sndBestName[ww[:,0], ww[:,1]],
    'tfidf': score[ww[:,0], ww[:,1]],
    'idf': idf[ww[:,0]],
    'qval': qvals[ww[:,0], ww[:,1]]
    })
    out = out.astype({'gene': 'str', 'cluster': 'str', 'secondBestClusterName': 'str'})
    return out        
