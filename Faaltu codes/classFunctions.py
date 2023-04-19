import numpy as np
import pandas as pd 
import statsmodels.api as sm
from scipy.stats import poisson

'''
This is an R function named estimateSoup which estimates the soup profile for a SoupChannel object.
The function takes in three arguments:

sc: a SoupChannel object
soupRange: a numeric vector of length 2 representing the range of droplet UMI counts to be considered in the soup estimation (default is c(0, 100))
keepDroplets: a logical value indicating whether or not to keep the droplet counts in the object after soup estimation (default is FALSE)

The function first checks if sc is a valid SoupChannel object using the is() function. Then, it estimates the soup profile by selecting droplets 
within the specified range (soupRange) and calculating the estimated soup expression for each gene using the rowSums() function. It saves the 
results in a data frame with row names from sc$tod and column names of "est" for estimated expression and "counts" for the total count of UMIs 
in each droplet. If keepDroplets is FALSE, the droplet count table (sc$tod) is removed to save space. The function returns the updated SoupChannel object.
'''
#tod -> each column is a droplet and each row is a gene
#toc -> just those columns of tod that contain a cell

def estimateSoup(sc, soupRange=[0, 100], keepDroplets=False):
    if not isinstance(sc, SoupChannel):
        raise ValueError("sc must be a SoupChannel object.")
    
    # Estimate the soup
    # w is a tuple of arrays, so we need to index the first element of the tuple to get the array we want
    w = np.where((sc.nDropUMIs > soupRange[0]) & (sc.nDropUMIs < soupRange[1]))[0]
    
    ''' Let the below columns be w in tod. Then,
    3 7 8   soup_est = 18/46, soup_counts = 18
    6 2 1   soup_est = 9/46, soup_counts = 9
    9 5 3   soup_est = 17/46, soup_counts = 17
    '''
    # np.sum(sc.tod[:, w], axis=1) computes the sum of the values in each row of the sc.tod array that correspond to the column indices specified in w. The 
    # resulting array will have the same number of rows as sc.tod. Dividing this sum by the total sum of the same columns, np.sum(sc.tod[:, w]), normalizes 
    # the result to be between 0 and 1. This gives a relative estimate of how much each row's "soup" quantity contributes to the overall "soup" in the selected columns.
    soup_est = np.sum(sc.tod[:, w], axis=1) / np.sum(sc.tod[:, w])
    soup_counts = np.sum(sc.tod[:, w], axis=1)

    soup_profile = pd.DataFrame({'est': soup_est, 'counts': soup_counts}, index=sc.tod.index)
    # index indicates the row names of the data frame corresponding to the genes in tod
    sc.soupProfile = soup_profile
    
    # Saves a lot of space if we can drop the droplets now we're done with them
    if not keepDroplets:
        sc.tod = None
    return sc

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

#tod -> each column is a droplet and each row is a gene
#toc -> just those columns of tod that contain a cell

def SoupChannel(tod, toc, metaData=None, calcSoupProfile=True, **kwargs):

    if metaData is not None and not all(np.sort(toc.columns) == np.sort(metaData.index)):
        raise ValueError("Rownames of metaData must match column names of table of counts.")
    # metadata ke row mai droplets hai, jinke andar cell hai
    
    if tod.shape[0] != toc.shape[0]:
        raise ValueError("The provided table of droplets (tod) and table of counts (toc) have different numbers of genes. Both tod and toc must have the same genes in the same order.")
    
    if not all(tod.index == toc.index):
        raise ValueError("Rownames of the table of droplets (tod) and table of counts (toc) differ. Both tod and toc must have the same genes in the same order.")
    
    # out is a dictionary
    out = {'tod': tod, 'toc': toc, **kwargs}
    
    metaData_df = pd.DataFrame(index=toc.columns, data={'nUMIs': toc.sum()})
    # metaData_df contains the number of UMIs associated with each droplet in toc. metaData_df ki index is equal to number of columns in toc
    out['metaData'] = metaData_df
        
    # This code is checking if the variable metaData is not None, which means that it has a value assigned to it. If metaData is not None, then it proceeds to 
    # execute the following two lines of code:
    #     metaData_df = metaData_df.drop(columns='nUMIs', errors='ignore')--> This line of code drops the column named 'nUMIs' from the metaData_df DataFrame. 
    #                 The errors='ignore' parameter is used to prevent an error from being raised if the 'nUMIs' column does not exist in the DataFrame.
    #     out['metaData'] = pd.concat([out['metaData'], metaData_df], axis=1, sort=False)--> This line of code concatenates the (new changed) metaData_df DataFrame with the 
    #                 'metaData' DataFrame (original dataframe) that is stored in the dictionary out. The resulting concatenated DataFrame is then assigned back to the 'metaData' key 
    #                 in the out dictionary.
    # In summary, this code updates the 'metaData' key in the out dictionary by dropping the 'nUMIs' column from the metaData_df DataFrame (if it exists), and 
    # concatenating the resulting DataFrame with the 'metaData' DataFrame that is already stored in the out dictionary.
    if metaData is not None:
        metaData_df = metaData_df.drop(columns='nUMIs', errors='ignore')
        out['metaData'] = pd.concat([out['metaData'], metaData_df], axis=1, sort=False)
    
    out['nDropUMIs'] = tod.sum()    # tod.sum() returns the sum of the values in each column of tod. The resulting array will have the same number of columns as tod.
    out['class'] = ['list', 'SoupChannel']  #indicates that the class of the object represented by the dictionary is a list of SoupChannel objects
    
    if calcSoupProfile:
        out = estimateSoup(out)     # out is sent as a SoupChannel object to the estimateSoup function
    ''' The dictionary out contains the following:
    tod : tod
    toc : toc
    metaData : metaData_df
    nDropUMIs : nDropUMIs
    class : ['list', 'SoupChannel']
    '''
    return out

def print_SoupChannel(x):
    print(f"Channel with {x['toc'].shape[0]} genes and {x['toc'].shape[1]} cells")