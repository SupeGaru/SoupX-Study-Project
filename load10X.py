# The R code defines a function load10X that loads unfiltered 10X data from each data-set and 
# identifies which droplets are cells using the cellranger defaults. It takes arguments like
#  dataDir which is the top level cellranger output directory (the directory that contains the raw_gene_bc_matrices folder),
#   cellIDs which are barcodes of droplets that contain cells and channelName which is the name of the channel to store. 
#   If cellIDs is NULL, the default cellranger set is used. Other arguments include readArgs which is a list of extra parameters
#    passed to Seurat::Read10X, includeFeatures which keeps only the types mentioned here and collapses to a single matrix,
#     verbose which is a logical variable for being verbose and ... which is for extra parameters passed to SoupChannel 
#     construction function.

# The function first works out which version of 10X data it's dealing with using the existence of specific directories. 
# It then loads the raw count data and identifies cell-only count data if cellIDs is not NULL. If it is NULL, 
# it works out which ones contain cells. The cluster annotation, fine-grained clusters, and tSNE projection are extracted 
# if available and returned as a SoupChannel object containing the count tables for the 10X dataset.

import os
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler

def load_10X(data_dir, cellIDs=None, channel_name=None, read_args=None, include_features=["Gene Expression"], verbose=True):
    #Work out which version we're dealing with
    isV3 = os.path.isdir(os.path.join(data_dir, "raw_feature_bc_matrix"))
    isV7 = os.path.isdir(os.path.join(data_dir, "analysis", "clustering", "gene_expression_graphclust"))
    isMulti = os.path.isdir(os.path.join(data_dir, "analysis", "clustering", "gex"))
    tgt = os.path.join(data_dir, "raw_feature_bc_matrix" if isV3 else "raw_gene_bc_matrices")
    #Add the reference genome for the non-V3 ones
    if not isV3:
        tgt = os.path.join(tgt, os.listdir(tgt)[0])
    if verbose:
        print("Loading raw count data")
    dat = pd.read_csv(tgt, sep="\t", header=None)
    if verbose:
        print("Loading cell-only count data")
    if cellIDs is not None:
        if not all(elem in dat.columns for elem in cellIDs):
            raise ValueError("Not all supplied cellIDs found in raw data.")
        dat_cells = dat[cellIDs]
    else:
        #Work out which ones contain cells
        tgt = os.path.join(data_dir, "filtered_feature_bc_matrix" if isV3 else "filtered_gene_bc_matrices")
        if not isV3:
            tgt = os.path.join(tgt, os.listdir(tgt)[0])
        dat_cells = pd.read_csv(tgt, sep="\t", header=None)
        #If it's a list of multiple types, have to decide what to include and collapse to one matrix.
        if isinstance(dat, dict):
            dat = pd.concat([dat[key] for key in includeFeatures], axis=1)
            dat_cells = pd.concat([dat_cells[key] for key in includeFeatures], axis=1)
    if verbose:
        print("Loading extra analysis data where available")
    #Get the cluster annotation if available
    mDat = None
    #What needs to be added to make V7 directory structure work
    v7Prefix = "gene_expression_" if isV7 else ""
    tgt = os.path.join(data_dir, "analysis", "clustering", "gex", "graphclust", "clusters.csv" if isMulti and not isV7 else os.path.join("analysis", "clustering", v7Prefix+"graphclust", "clusters.csv"))
    if os.path.isfile(tgt):
        clusters = pd.read_csv(tgt)
        mDat = pd.DataFrame({"clusters": clusters["Cluster"]}, index=clusters["Barcode"])
    #Add fine grained clusters too if present
    tgt = os.path.join(data_dir, "analysis", "clustering", "gex", "kmeans_10_clusters", "clusters.csv" if isMulti and not isV7 else os.path.join("analysis", "clustering", v7Prefix+"kmeans_10_clusters", "clusters.csv"))
    if os.path.isfile(tgt):
        clusters = pd.read_csv(tgt)
        mDat["clustersFine"] = clusters["Cluster"]
         # Get tSNE if available and point to it
    if isMulti and not isV7:
        tgt = os.path.join(dataDir, 'analysis', 'dimensionality_reduction', 'gex', 'tsne_projection.csv')
    else:
        tgt = os.path.join(dataDir, 'analysis', 'tsne', f"{v7Prefix}2_components", 'projection.csv')

    if os.path.exists(tgt):
        tsne = pd.read_csv(tgt)
        if mDat is None:
            mDat = pd.DataFrame({'tSNE1': tsne['TSNE.1'], 'tSNE2': tsne['TSNE.2']}, index=tsne['Barcode'])
        else:
            mDat['tSNE1'] = tsne['TSNE.1'][tsne['Barcode'].isin(mDat.index)]
            mDat['tSNE2'] = tsne['TSNE.2'][tsne['Barcode'].isin(mDat.index)]
        DR = ['tSNE1', 'tSNE2']
    else:
        DR = None

    # Ensure rownames of metadata match column names of counts
    if mDat is not None and any(mDat.index != datCells.columns):
        mDat.index = mDat.index.str.replace('-1$', '')
        if any(mDat.index != datCells.columns):
            raise ValueError("Error matching meta-data to cell names.")

    # Get a name for the channel
    if channelName is None:
        channelName = dataDir if names(dataDir) is None else names(dataDir)

    return SoupChannel(tod=dat,
                       toc=datCells,
                       metaData=mDat,
                       channelName=channelName,
                       dataDir=dataDir,
                       dataType='10X',
                       isV3=not isV7,
                       DR=DR)