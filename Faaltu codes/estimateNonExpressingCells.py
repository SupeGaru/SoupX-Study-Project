'''
The R function estimateNonExpressingCells takes in a SoupChannel object sc, a list of sets of non-expressed genes nonExpressedGeneList, optional cluster assignments
clusters, maximum contamination fraction maximumContamination, and false discovery rate FDR. It returns a matrix of binary values indicating whether each cell is 
non-expressing or not.

If clusters is not provided or is set to NULL, the function will use individual cells as their own clusters. The function checks that all column names in the table 
of counts appear in the clusters vector. If nonExpressedGeneList is not a list, the function will throw an error.

The function then selects which genes to use for each cluster by finding unique genes in the nonExpressedGeneList and then selecting cells with the count data of 
these genes. The expected number of counts from contamination is then calculated, and the cells where a gene is definitely expressed are identified. This is done by 
applying the ppois function with the lower.tail=FALSE argument and then adjusting the p-values using the Benjamini-Hochberg method. The resulting matrix is returned.

If no non-expressing cells are identified, the function returns a warning suggesting that clusters be set to FALSE, or maximumContamination and/or FDR be increased. 
If fewer than 100 non-expressing cells are identified, the function returns a warning that the estimation of contamination fraction may be inaccurate.
'''

def estimate_non_expressing_cells(sc, non_expressed_gene_list, clusters=None, maximum_contamination=1.0, FDR=0.05):
    import numpy as np
    import pandas as pd
    import classFunctions as cf
    import statsmodels.api as sm
    from scipy.stats import poisson

    if not isinstance(sc, cf.SoupChannel):
        raise ValueError("sc is not a valid SoupChannel object")

    # Get clusters if they exist, if they don't, set to individual cells
    if clusters is None:
        if "clusters" in sc.metaData.columns:
            clusters = sc.metaData["clusters"].astype(str).to_dict()

    # Using each cell as its own cluster
    if clusters is None or (len(clusters) == 1 and list(clusters.values())[0] == "False"):
        print("No clusters found or supplied, using every cell as its own cluster.")
        clusters = dict(zip(sc.metaData.index, sc.metaData.index))

    # Check we have coverage of everything
    if not all(sc.toc.columns.isin(clusters.keys())):
        raise ValueError("Invalid cluster specification. clusters must be a named dict with all column names in the table of counts appearing.")

    # Convert gene list to genuine list if vector
    if not isinstance(non_expressed_gene_list, list):
        raise ValueError("non_expressed_gene_list must be a list of sets of genes. e.g. [{'HB': ['HBB', 'HBA2']}]")

    # Now work out which clusters to use which genes on
    tgt_gns = np.unique([gene for gene_set in non_expressed_gene_list for gene in gene_set.values()])
    dat = sc.toc.loc[tgt_gns]
    cnts = pd.concat([dat[e].sum(axis=1) for e in non_expressed_gene_list], axis=1, keys=non_expressed_gene_list.keys())
    
    # Work out how many counts we'd expect if the cell were maximally contaminated and all expression came from the contamination
    exp = np.outer(sc.soupProfile.loc[tgt_gns, "est"], sc.metaData["nUMIs"] * maximum_contamination)
    exp = pd.DataFrame(exp, index=tgt_gns, columns=sc.metaData.index)
    exp = pd.concat([exp[e].sum(axis=1) for e in non_expressed_gene_list], axis=1, keys=non_expressed_gene_list.keys())

    # Identify those cells where a gene is definitely expressed
    s = pd.Series(clusters)
    clust_exp = pd.DataFrame(poisson.pmf(cnts - 1, exp), index=tgt_gns, columns=cnts.columns)
    clust_exp = clust_exp.apply(lambda x: x.apply(lambda y: min(y)), axis=1)
    clust_exp = clust_exp.apply(lambda x: x.to_numpy(), axis=1, result_type="broadcast")
    clust_exp = pd.concat([clust_exp[s[s == c].index] for c in s.unique()], axis=1, keys=s.unique())
    clust_exp = clust_exp.apply(lambda x: x >= FDR)

    # Expand it out into a full cell matrix
    clust_exp = clust_exp.loc[s.index]
    clust_exp.index = s.index

    # Check that we actually got some
    if clust_exp.sum().sum() == 0:
        print("No non-expressing cells identified. Consider setting clusters=FALSE, increasing maximumContamination and/or FDR")
    
    # A small number found
    if 0 < clust_exp.sum().sum() < 100:
        print("Fewer than 100 non-expressing cells identified.  The estimation of the contamination fraction may be inaccurate.  Consider setting clusters=FALSE, increasing maximumContamination and/or FDR")

    # Return the results
    return clust_exp