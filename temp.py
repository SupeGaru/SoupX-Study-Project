import numpy as np
import pandas as pd 
import statsmodels.api as sm
from scipy.stats import poisson
import matplotlib.pyplot as pl

tod = pd.read_csv("data.csv", index_col=0)
tod.drop('Barcode', axis=1, inplace=True)
empty_drops = tod.loc[tod['sum']<20]
toc = tod.loc[tod['sum']>20]

def auxiliary():
    geneColumns = ['swH3N2_PB2_comm_raw', 'swH1N1_PB2_comm_raw', 'swH3N2_PB1_comm2_raw', 'swH1N1_PB1_comm2_raw',
                'swH3N2_PA_comm_raw', 'swH1N1_PA_comm_raw', 'swH3N2_HAa_raw', 'swH1N1_HAa_raw', 
                'swH3N2_NP_comm_raw', 'swH1N1_NP_comm_raw', 'swH3N2_NAa_raw', 'swH1N1_NAb_raw',
                'swH3N2_M_comm_raw', 'swH1N1_M_comm_raw', 'swH3N2_NS_comm_raw', 'swH1N1_NS_comm_raw']
    detected = toc[['subsets_H1N1_sum', 'subsets_H3N2_sum']]
    h1n1_columns = toc[['swH1N1_PB2_comm_raw', 'swH1N1_PB1_comm2_raw', 'swH1N1_PA_comm_raw', 'swH1N1_HAa_raw',
                    'swH1N1_NP_comm_raw', 'swH1N1_NAb_raw', 'swH1N1_M_comm_raw', 'swH1N1_NS_comm_raw']]
    h3n2_columns = toc[['swH3N2_PB2_comm_raw', 'swH3N2_PB1_comm2_raw', 'swH3N2_PA_comm_raw', 'swH3N2_HAa_raw',
                    'swH3N2_NP_comm_raw', 'swH3N2_NAa_raw', 'swH3N2_M_comm_raw', 'swH3N2_NS_comm_raw']]

    h1n1_column_names = ['swH1N1_PB2_comm_raw', 'swH1N1_PB1_comm2_raw', 'swH1N1_PA_comm_raw', 'swH1N1_HAa_raw',
                    'swH1N1_NP_comm_raw', 'swH1N1_NAb_raw', 'swH1N1_M_comm_raw', 'swH1N1_NS_comm_raw']
    h3n2_column_names = ['swH3N2_PB2_comm_raw', 'swH3N2_PB1_comm2_raw', 'swH3N2_PA_comm_raw', 'swH3N2_HAa_raw',
                    'swH3N2_NP_comm_raw', 'swH3N2_NAa_raw', 'swH3N2_M_comm_raw', 'swH3N2_NS_comm_raw']
    return geneColumns, detected, h1n1_columns, h3n2_columns, h1n1_column_names, h3n2_column_names

def create_bg_gene_empty():
    geneColumns = auxiliary()[0]
    #sum of an individual gene in all empty drops , size = number of genes
    geneSum_empty = np.sum(empty_drops[geneColumns], axis=0)    # geneSum_empty is list
    #sum of all data
    totalSum_empty = np.sum(geneSum_empty)                      # totalSum_empty is a number
    #fraction of background expression of all genes ,size = number of genes
    bg_gene_empty = geneSum_empty/totalSum_empty                # bg_gene_empty (for empty droplets) is list
    return bg_gene_empty
    # equation 1 in paper done, bg_gene is bg in formula

#sum of individual genes in all toc drops, size = number of empty drops
# geneSum_toc = np.sum(toc[geneColumns], axis=0)              # geneSum_toc is list

# found the set of genes/cells for which we can assume that there is no cell endogenous expression (wahan nahi hone chahiye tha) (equation 4 ke upar)

def identify_doublets():
    geneColumns = auxiliary()[0]
    detected = auxiliary()[1]
    global toc

    for i in range(len(detected)):
        h1n1 = detected['subsets_H1N1_sum'][i]
        h3n2 = detected['subsets_H3N2_sum'][i]

        if(h3n2 > h1n1 and h3n2/(h3n2+h1n1)<0.8):
            for j in range(len(geneColumns)):
                toc[geneColumns[j]][i] = 0.0001
        elif(h1n1 > h3n2 and h1n1/(h3n2+h1n1)<0.8):
            for j in range(len(geneColumns)):
                toc[geneColumns[j]][i] = 0.0001
    

def create_rho_c():
    detected = auxiliary()[1]
    h1n1_columns = auxiliary()[2]
    h3n2_columns = auxiliary()[3]
    h1n1_column_names = auxiliary()[4]
    h3n2_column_names = auxiliary()[5]
    bg_gene_empty = create_bg_gene_empty()

    rho_c = []

    for i in range(len(detected)):
        h1n1 = detected['subsets_H1N1_sum'][i]
        h3n2 = detected['subsets_H3N2_sum'][i]

        ngc_sum = 0
        nc = 0
        bgc_sum = 0

        if(h1n1 < h3n2):
            ngc_sum = h1n1
            nc = h3n2 + h1n1
            if(h1n1>0):
                for k in range(len(h1n1_column_names)):
                    if(h1n1_columns[h1n1_column_names[k]][i])>0:
                        bgc_sum += bg_gene_empty[2*k+1]
                rho_c.append(ngc_sum/(nc*bgc_sum))
            else:
                rho_c.append(0)

        elif(h1n1 > h3n2):
            ngc_sum = h3n2
            nc = h3n2 + h1n1
            if(h3n2>0):
                for k in range(len(h3n2_column_names)):
                    if(h3n2_columns[h3n2_column_names[k]][i])>0:
                        bgc_sum += bg_gene_empty[2*k]
                rho_c.append(ngc_sum/(nc*bgc_sum))
            else:
                rho_c.append(0)
        
        else:
            rho_c.append(0)

    return rho_c

def create_ogc():
    rho_c = create_rho_c()
    bg_gene_empty = create_bg_gene_empty()
    nc_total = toc['sum']
    ogc_rows, ogc_cols = (len(rho_c), len(bg_gene_empty))
    ogc = [[0 for i in range(ogc_cols)] for j in range(ogc_rows)]

    for i in range(len(rho_c)):
        for j in range(len(bg_gene_empty)):
            ogc[i][j] = rho_c[i]*bg_gene_empty[j]*nc_total[i]
    # print(ogc)
    return ogc

def create_mgc(toc):
    bg_gene_empty = create_bg_gene_empty()
    ogc = create_ogc()
    geneColumns = auxiliary()[0]
    ogc_rows, ogc_cols = (len(ogc), len(ogc[0]))
    ngc = toc[geneColumns]
    mgc = [[0 for i in range(ogc_cols)] for j in range(ogc_rows)]
    for i in range(len(ngc)):
        for j in range(len(ngc.columns)):
            mgc[i][j] = ngc.iloc[i][j] - ogc[i][j]
            if(mgc[i][j]<0):
                mgc[i][j] = 0
    return mgc

def write_mgc():
    mgc = create_mgc(toc)
    mgc = pd.DataFrame(mgc)
    mgc.to_csv('output.csv', index=False, header=False)

identify_doublets()
write_mgc()