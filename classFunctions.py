'''
The SoupChannel function takes in two tables of data: tod, which is a table of droplet barcodes and their associated gene names, and toc, which is
a table of gene counts for each droplet.
The function checks that the dimensions and row/column names of tod and toc are compatible, and then creates a list object out that contains both
tables and any additional arguments that are passed to the function using the ... syntax.

If metaData is provided, the function checks that the row names of metaData match the column names of toc, and then creates a data frame 
metaData_df that contains the number of unique molecular identifiers (UMIs) for each barcode in toc. If metaData is not provided, metaData_df is 
simply the nUMIs column of toc. metaData is then merged into metaData_df, if it is not NULL.

The function calculates the number of UMIs associated with each droplet barcode in tod and adds it to the output list as nDropUMIs.
If calcSoupProfile is TRUE, the function calls the estimateSoup function on the output list object.
The print.SoupChannel function simply prints a summary of the SoupChannel object, including the number of genes and cells in the toc table.

Overall, the SoupChannel function is used to create a SoupChannel object, which contains information about droplet barcodes, gene counts, and any
additional metadata associated with the data. The estimateSoup function is used to estimate the "soup" of background RNA that contaminates droplets
in single-cell RNA sequencing experiments.
'''

def SoupChannel(tod, toc, metaData=None, calcSoupProfile=True, **kwargs):
    import estimateSoup
    import numpy as np
    import pandas as pd
    import statsmodels.api as sm
    from scipy.stats import poisson
    
    if metaData is not None and not all(np.sort(toc.columns) == np.sort(metaData.index)):
        raise ValueError("Rownames of metaData must match column names of table of counts.")
    
    if tod.shape[0] != toc.shape[0]:
        raise ValueError("The provided table of droplets (tod) and table of counts (toc) have different numbers of genes. Both tod and toc must have the same genes in the same order.")
    
    if not all(tod.index == toc.index):
        raise ValueError("Rownames of the table of droplets (tod) and table of counts (toc) differ. Both tod and toc must have the same genes in the same order.")
    
    out = {'tod': tod, 'toc': toc, **kwargs}
    
    metaData_df = pd.DataFrame(index=toc.columns, data={'nUMIs': toc.sum()})
    out['metaData'] = metaData_df
    
    if metaData is not None:
        metaData_df = metaData_df.drop(columns='nUMIs', errors='ignore')
        out['metaData'] = pd.concat([out['metaData'], metaData_df], axis=1, sort=False)
    
    out['nDropUMIs'] = tod.sum()
    out['class'] = ['list', 'SoupChannel']
    
    if calcSoupProfile:
        out = estimateSoup(out)
    
    return out


def print_SoupChannel(x):
    print(f"Channel with {x['toc'].shape[0]} genes and {x['toc'].shape[1]} cells")

# Note: You will need to define the estimateSoup function before using the SoupChannel function in Python.