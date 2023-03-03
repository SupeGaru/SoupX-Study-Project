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

def estimateSoup(sc, soupRange=[0, 100], keepDroplets=False):
    import numpy as np
    import pandas as pd
    import classFunctions as cf
    import statsmodels.api as sm
    #from scipy.stats import poisson

    if not isinstance(sc, cf.SoupChannel):
        raise ValueError("sc must be a SoupChannel object.")
    
    # Estimate the soup
    w = np.where((sc.nDropUMIs > soupRange[0]) & (sc.nDropUMIs < soupRange[1]))[0]
    soup_est = np.sum(sc.tod[:, w], axis=1) / np.sum(sc.tod[:, w])
    soup_counts = np.sum(sc.tod[:, w], axis=1)
    soup_profile = pd.DataFrame({'est': soup_est, 'counts': soup_counts}, index=sc.tod.index)
    sc.soupProfile = soup_profile
    
    # Saves a lot of space if we can drop the droplets now we're done with them
    if not keepDroplets:
        sc.tod = None
    return sc

# Note that this code assumes that the SoupChannel object has the following attributes: nDropUMIs, tod, and soupProfile.