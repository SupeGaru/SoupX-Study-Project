# Python equivalent code for creating PBMC_metaData object
import numpy as np
import soupX
import pandas as pd
import scanpy as sc
from sklearn.preprocessing import StandardScaler

# Load PBMC 4k data
adata = sc.read_10x_mtx("path/to/10x/directory")

# Preprocessing steps
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
sc.pp.pca(adata, n_comps=30, use_highly_variable=True)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=1)

# Create PBMC_metaData object
PBMC_metaData = pd.DataFrame(adata.obsm["X_umap"], columns=["RD1", "RD2"])
PBMC_metaData["Cluster"] = adata.obs["leiden"]
PBMC_metaData["Annotation"] = adata.obs["leiden"].replace(
    {
        "0": "MNP",
        "1": "T_CD4",
        "2": "T_CD4",
        "3": "T_CD8",
        "4": "B",
        "5": "T_CD8",
        "6": "NK",
        "7": "B",
        "8": "NK",
        "9": "MNP",
        "10": "MNP",
        "11": "?",
    }
)

# Python equivalent code for creating PBMC_sc object

# Load PBMC 4k data
data_dir = "path/to/10x/directory"
PBMC_sc = soupX.load10X(data_dir, calcSoupProfile=False)
PBMC_sc = soupX.SoupChannel(
    PBMC_sc.tod, PBMC_sc.toc[:, list(range(0, 2170, 2))])

# Python equivalent code for creating scToy object

scToy = soupX.SoupChannel(
    np.random.rand(62, 226),
    cluster=np.random.choice(range(5), size=62),
    rd1=np.random.normal(size=62),
    rd2=np.random.normal(size=62),
    contamination_estimate=0.1,
)
