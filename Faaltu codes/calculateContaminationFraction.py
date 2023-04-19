"""
The given R code defines a function called calculateContaminationFraction that takes four arguments:

sc: a SoupChannel object
nonExpressedGeneList: a list of sets of genes
useToEst: a logical matrix indicating which cells should be used for estimation
verbose: a logical indicating whether to print progress messages (default: TRUE)
forceAccept: a logical indicating whether to force acceptance of the estimated contamination fraction (default: FALSE)
The function first checks that sc is of the correct type, and that nonExpressedGeneList is a list of sets of genes. It also checks that at least one cell is specified as acceptable for estimation.

The function then constructs a data frame based on the non-expressed gene sets, and fits a Poisson GLM with a log-link to estimate the global contamination fraction. The function then adds the estimated contamination fraction and its confidence interval to the meta-data of sc.
Finally, the function returns the modified sc object.
The Python code is a translation of the R code, using pandas and statsmodels instead of data frames and glm in R. The only differences are in the syntax and package functions used.
"""

def calculate_contamination_fraction(sc, non_expressed_gene_list, use_to_est, verbose=True, force_accept=False):
    import numpy as np
    import pandas as pd
    import setProperties as sp
    import classFunctions as cf
    import statsmodels.api as sm
    from scipy.stats import poisson
    
    if not isinstance(sc, cf.SoupChannel):
        raise ValueError("sc must be a SoupChannel object")
    
    if not isinstance(non_expressed_gene_list, list):
        raise ValueError("nonExpressedGeneList must be a list of sets of genes. e.g. [{'HB': ['HBB', 'HBA2']}]")
    
    if sum(use_to_est) == 0:
        raise ValueError("No cells specified as acceptable for estimation. useToEst must not be all False")
    
    df = []
    for i, gene_set in enumerate(non_expressed_gene_list):
        tgts = gene_set.values()
        sFrac = np.sum(sc.soupProfile.loc[tgts, 'est'])
        w = use_to_est[:, i].index[use_to_est[:, i]]
        if len(w) > 0:
            cnts = sc.toc.loc[tgts, w]
            df.append(pd.DataFrame({'cells': cnts.columns, 'geneSet': i, 'soupFrac': sFrac,
                                    'counts': cnts.sum(axis=0)}, index=range(cnts.shape[1])))
    if len(df) > 0:
        df = pd.concat(df)
        df['nUMIs'] = sc.metaData.loc[df['cells'], 'nUMIs'].values
        df['expSoupCnts'] = df['nUMIs'] * df['soupFrac']
        df['counts'] = df['counts'].astype(int)
        # Make cells a categorical variable, but preserve ordering in sc. 
        # This ensures that the STAN index is the index in sc.metaData.
        df['cells'] = pd.Categorical(df['cells'], categories=sc.metaData.index.tolist())
        fit = sm.GLM(df['counts'], sm.add_constant(np.zeros(df.shape[0])), family=poisson(),
                     offset=np.log(df['expSoupCnts'])).fit()
        # Finally, add it to the meta-data
        sc = sp.set_contamination_fraction(sc, np.exp(fit.params[0]), force_accept=force_accept)
        ci = fit.conf_int(alpha=0.05)
        sc.metaData['rhoLow'] = np.exp(ci.iloc[0, 0])
        sc.metaData['rhoHigh'] = np.exp(ci.iloc[0, 1])
        if verbose:
            print(f"Estimated global contamination fraction of {100 * np.exp(fit.params[0]):.2f}%")
    return sc

"""
Note that in the above Python code, I had to make a few assumptions about the definitions of SoupChannel, set_contamination_fraction(), 
and some of the data structures used in the R code. You may need to adjust the code slightly to match your specific use case.
"""