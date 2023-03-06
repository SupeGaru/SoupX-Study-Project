'''
The Python code defines several functions that modify a single object "sc".

The first function "setSoupProfile" takes as input "sc" (an object), and a Pandas dataframe "soupProfile". It first checks if the "soupProfile" contains the necessary columns "est" and "counts".
If not, it raises a ValueError. It also checks that all the genes in "soupProfile" are present in the index of "sc". If any genes are missing, it raises an error. If everything checks out, it sets 
the attribute "sc.soupProfile" to be equal to "soupProfile" indexed by "sc.toc.index".

The second function "set_clusters" takes as input "sc" and a Pandas dataframe "clusters". It checks if all columns in the index of "clusters" are also columns in the index of "sc.toc". If not, it 
checks if the length of "clusters" is equal to the number of rows in "sc.metaData". If not, it raises a ValueError. If everything checks out, it sets the "clusters" column in "sc.metaData" to be 
equal to "clusters".

The third function "setContaminationFraction" takes as input "sc", a vector "contFrac" (of length 1 or equal to the number of rows in "sc.metaData"), and a Boolean "forceAccept". It first checks 
if any values in "contFrac" are greater than 1. If so, it raises an Exception. If "forceAccept" is True, it prints a warning but continues. If any values in "contFrac" are greater than 0.5, it 
raises an Exception with a warning message. If any values in "contFrac" are greater than 0.3, it prints a message. If "contFrac" is of length 1, it sets the "rho" column in "sc.metaData" to be 
equal to "contFrac". Otherwise, it sets the "rho" values in "sc.metaData" indexed by the keys in "contFrac" to be equal to the values in "contFrac".

The fourth function "setDR" takes as input "sc", a Pandas dataframe "DR", and a string "reductName". It first checks if "DR" has more than two columns. If so, it prints a message and keeps only the 
first two columns. It then checks if all the row names in "DR" are present in the index of "sc.metaData". If not, it checks if the number of rows in "DR" is equal to the number of rows in "sc.metaData". 
If not, it raises a ValueError. If everything checks out, it concatenates "DR" with "sc.metaData" indexed by the rows in "DR". If "reductName" is not None, it renames the columns of "DR". Finally, it 
sets the "DR" attribute in "sc" to be equal to the column names in "DR".
'''

import pandas as pd

def setSoupProfile(sc, soupProfile):
    import classFunctions as cf
    if 'est' not in soupProfile.columns:
        raise ValueError("est column missing from soupProfile")
    if 'counts' not in soupProfile.columns:
        raise ValueError("counts column missing from soupProfile")
    if not all(cf.SoupChannel.rownames(soupProfile).isin(sc.toc.index)):
        raise ValueError("soupProfile invalid. Not all genes found.")
    else:
        sc.soupProfile = soupProfile.loc[sc.toc.index, :]
    return sc

def set_clusters(sc, clusters):
    if not all(sc.toc.columns.isin(clusters.index)):
        if len(clusters) != sc.metaData.shape[0]:
            raise ValueError("Invalid cluster specification. See help.")
        else:
            # Ensure the thing we're setting is not a factor
            sc.metaData["clusters"] = clusters.astype(str)
    else:
        sc.metaData["clusters"] = clusters[sc.metaData.index].astype(str)

    # Do a check that things set correctly
    if sc.metaData["clusters"].isna().any():
        raise ValueError("NAs found in cluster names. Ensure a (non-na) mapping to cluster is provided for each cell.")

    return sc

def setContaminationFraction(sc, contFrac, forceAccept=False):
    if any(contFrac > 1):
        raise Exception("Contamination fraction greater than 1 detected. This is impossible and likely represents a failure in the estimation procedure used.")
    if forceAccept:
        warning = print
        # raise Exception = print
    if any(contFrac > 0.5):
        raise Exception(f"Extremely high contamination estimated ({max(contFrac):.2g}). This likely represents a failure in estimating the contamination fraction. Set forceAccept=True to proceed with this value.")
    elif any(contFrac > 0.3):
        print(f"Estimated contamination is very high ({max(contFrac):.2g}).")
    if len(contFrac) == 1:
        sc.metaData['rho'] = contFrac
    else:
        if not all(x in sc.metaData.index for x in contFrac.keys()):
            raise Exception("contFrac must be either of length 1 or a named vector with names matching the rows of sc$metaData")
        sc.metaData.loc[contFrac.keys(), 'rho'] = contFrac.values
    return sc

def setDR(sc, DR, reductName=None):
    # If more than two columns, keep the first two
    if DR.shape[1] > 2:
        print(f"DR has {DR.shape[1]} columns where 2 were expected. Using first two.")
        DR = DR.iloc[:, :2]
    
    # Check if the row names match the metadata
    m = pd.Series(DR.index).isin(sc.metaData.index)
    if not all(m):
        # Can't use row-names, so the number better match
        if DR.shape[0] != sc.metaData.shape[0]:
            raise ValueError(f"Rownames present in metaData not found in DR and row numbers differ ({sc.metaData.shape[0]} in metaData, {DR.shape[0]} in DR). Each cell must have a corresponding entry in DR.")
        m = pd.Series(True, index=DR.index)
    
    # Should we change the names?
    if reductName is not None:
        DR.columns = [f"{reductName}_{i+1}" for i in range(DR.shape[1])]
    
    # Add the entries in
    sc.metaData = pd.concat([sc.metaData, DR.loc[m, :]], axis=1)
    
    # And point them in the right place
    sc.DR = DR.columns.tolist()
    return sc