import numpy as np
import pandas as pd 
import statsmodels.api as sm
from scipy.stats import poisson
import matplotlib.pyplot as plt

tod = pd.read_csv("data.csv", index_col=0)
# tod.drop('Barcode', axis=1, inplace=True)
empty_drops = tod.loc[tod['sum']<20]
toc = tod.loc[tod['sum']>20]
barcode = toc['Barcode']
toc.drop('Barcode', axis=1, inplace=True)
xInitialSum = []
yInitialSum = []
errorInitial = 0
errorFinal = 0

# print("Size of toc is: ", len(toc), "\n")

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

#sum of individual genes in all toc drops, size = number of empty drops
# geneSum_toc = np.sum(toc[geneColumns], axis=0)              # geneSum_toc is list

# found the set of genes/cells for which we can assume that there is no cell endogenous expression (wahan nahi hone chahiye tha) (equation 4 ke upar)

def create_rho_c():
    detected = auxiliary()[1]
    h1n1_columns = auxiliary()[2]
    h3n2_columns = auxiliary()[3]
    h1n1_column_names = auxiliary()[4]
    h3n2_column_names = auxiliary()[5]
    bg_gene_empty = create_bg_gene_empty()

    global xInitialSum, yInitialSum, errorInitial, errorFinal

    rho_c = []

    for i in range(len(detected)):
        h1n1 = detected['subsets_H1N1_sum'][i]
        h3n2 = detected['subsets_H3N2_sum'][i]

        xInitialSum.append(h1n1)
        yInitialSum.append(h3n2)

        ngc_sum = 0
        nc = 0
        bgc_sum = 0

        if(h1n1 < h3n2):
            errorInitial += h1n1
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
            errorInitial += h3n2
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
    old_mgc = create_mgc(toc)
    new_mgc = []
    headers = ['Barcode', 'swH3N2_PB2_comm_raw', 'swH1N1_PB2_comm_raw', 'swH3N2_PB1_comm2_raw', 'swH1N1_PB1_comm2_raw',
                'swH3N2_PA_comm_raw', 'swH1N1_PA_comm_raw', 'swH3N2_HAa_raw', 'swH1N1_HAa_raw', 
                'swH3N2_NP_comm_raw', 'swH1N1_NP_comm_raw', 'swH3N2_NAa_raw', 'swH1N1_NAb_raw',
                'swH3N2_M_comm_raw', 'swH1N1_M_comm_raw', 'swH3N2_NS_comm_raw', 'swH1N1_NS_comm_raw']
    new_mgc.append(headers)
    for i in range(len(old_mgc)):
        temp = []
        temp.append(barcode[i])
        # new_mgc.append(barcode[i])
        for j in range(len(old_mgc[0])):
            temp.append(old_mgc[i][j])
        new_mgc.append(temp)
    mgc = pd.DataFrame(new_mgc)
    mgc.to_csv('output.csv', index=False, header=False)

    return old_mgc

def SumPlotGraph():
    mgc = write_mgc()
    geneColumns = auxiliary()[0]
    xFinalSum = []
    yFinalSum = []

    for i in range(len(mgc)):
        h1n1 = 0
        h3n2 = 0
        for j in range (len(geneColumns)):
            if j%2==0 :
                h3n2 += mgc[i][j]
            else:
                h1n1 += mgc[i][j]
        
        xFinalSum.append(h1n1)
        yFinalSum.append(h3n2)

    plt.subplot(1,2,1)
    plt.plot(xInitialSum, yInitialSum, 'o', color = 'blue')
    plt.title('Initial nUMIs')
    plt.xlabel('H1N1')
    plt.ylabel('H3N2')

    plt.subplot(1,2,2)
    plt.plot(xFinalSum, yFinalSum, 'o', color = 'red')
    plt.title('Final nUMIs')
    plt.xlabel('H1N1')
    plt.ylabel('H3N2')
    plt.show()

def individualInitialPlotGraph():
    # here mgc is toc, used just coz code copy kiya hai neeche se
    mgc = toc[['swH3N2_PB2_comm_raw', 'swH1N1_PB2_comm_raw', 'swH3N2_PB1_comm2_raw', 'swH1N1_PB1_comm2_raw',
                'swH3N2_PA_comm_raw', 'swH1N1_PA_comm_raw', 'swH3N2_HAa_raw', 'swH1N1_HAa_raw', 
                'swH3N2_NP_comm_raw', 'swH1N1_NP_comm_raw', 'swH3N2_NAa_raw', 'swH1N1_NAb_raw',
                'swH3N2_M_comm_raw', 'swH1N1_M_comm_raw', 'swH3N2_NS_comm_raw', 'swH1N1_NS_comm_raw']]

    # x is H1N1, y is H3N2
    xPB2 = []
    yPB2 = []
    xPB1 = []
    yPB1 = []
    xPA = []
    yPA = []
    xHAa = []
    yHAa = []
    xNP = []
    yNP = []
    yNAa = []
    xNAb = []
    xM = []
    yM = []
    xNS = []
    yNS = []

    for i in range (len(mgc)):
        yPB2.append(mgc.iloc[i][0])
        xPB2.append(mgc.iloc[i][1])
        yPB1.append(mgc.iloc[i][2])
        xPB1.append(mgc.iloc[i][3])
        yPA.append(mgc.iloc[i][4])
        xPA.append(mgc.iloc[i][5])
        yHAa.append(mgc.iloc[i][6])
        xHAa.append(mgc.iloc[i][7])
        yNP.append(mgc.iloc[i][8])
        xNP.append(mgc.iloc[i][9])
        yNAa.append(mgc.iloc[i][10])
        xNAb.append(mgc.iloc[i][11])
        yM.append(mgc.iloc[i][12])
        xM.append(mgc.iloc[i][13])
        yNS.append(mgc.iloc[i][14])
        xNS.append(mgc.iloc[i][15])

    plt.subplot(4,2,1)
    plt.plot(xPB2, yPB2, 'o', color = 'blue')
    plt.title('PB2')
    plt.xlabel('H1N1')
    plt.ylabel('H3N2')

    plt.subplot(4,2,2)
    plt.plot(xPB1, yPB1, 'o', color = 'blue')
    plt.title('PB1')
    plt.xlabel('H1N1')
    plt.ylabel('H3N2')

    plt.subplot(4,2,3)
    plt.plot(xPA, yPA, 'o', color = 'blue')
    plt.title('PA')
    plt.xlabel('H1N1')
    plt.ylabel('H3N2')

    plt.subplot(4,2,4)
    plt.plot(xHAa, yHAa, 'o', color = 'blue')
    plt.title('HAa')
    plt.xlabel('H1N1')
    plt.ylabel('H3N2')

    plt.subplot(4,2,5)
    plt.plot(xNP, yNP, 'o', color = 'blue')
    plt.title('NP')
    plt.xlabel('H1N1')
    plt.ylabel('H3N2')

    plt.subplot(4,2,6)
    plt.plot(xNAb, yNAa, 'o', color = 'blue')
    plt.title('NAa/b')
    plt.xlabel('H1N1')
    plt.ylabel('H3N2')

    plt.subplot(4,2,7)
    plt.plot(xM, yM, 'o', color = 'blue')
    plt.title('M')
    plt.xlabel('H1N1')
    plt.ylabel('H3N2')

    plt.subplot(4,2,8)
    plt.plot(xNS, yNS, 'o', color = 'blue')
    plt.title('NS')
    plt.xlabel('H1N1')
    plt.ylabel('H3N2')

    plt.show()

def individualFinalPlotGraph():
    mgc = write_mgc()
    # x is H1N1, y is H3N2
    xPB2 = []
    yPB2 = []
    xPB1 = []
    yPB1 = []
    xPA = []
    yPA = []
    xHAa = []
    yHAa = []
    xNP = []
    yNP = []
    yNAa = []
    xNAb = []
    xM = []
    yM = []
    xNS = []
    yNS = []


    for i in range(len(mgc)):
        yPB2.append(mgc[i][0])
        xPB2.append(mgc[i][1])
        yPB1.append(mgc[i][2])
        xPB1.append(mgc[i][3])
        yPA.append(mgc[i][4])
        xPA.append(mgc[i][5])
        yHAa.append(mgc[i][6])
        xHAa.append(mgc[i][7])
        yNP.append(mgc[i][8])
        xNP.append(mgc[i][9])
        yNAa.append(mgc[i][10])
        xNAb.append(mgc[i][11])
        yM.append(mgc[i][12])
        xM.append(mgc[i][13])
        yNS.append(mgc[i][14])
        xNS.append(mgc[i][15])

    plt.subplot(4,2,1)
    plt.plot(xPB2, yPB2, 'o', color = 'blue')
    plt.title('PB2')
    plt.xlabel('H1N1')
    plt.ylabel('H3N2')

    plt.subplot(4,2,2)
    plt.plot(xPB1, yPB1, 'o', color = 'blue')
    plt.title('PB1')
    plt.xlabel('H1N1')
    plt.ylabel('H3N2')

    plt.subplot(4,2,3)
    plt.plot(xPA, yPA, 'o', color = 'blue')
    plt.title('PA')
    plt.xlabel('H1N1')
    plt.ylabel('H3N2')

    plt.subplot(4,2,4)
    plt.plot(xHAa, yHAa, 'o', color = 'blue')
    plt.title('HAa')
    plt.xlabel('H1N1')
    plt.ylabel('H3N2')

    plt.subplot(4,2,5)
    plt.plot(xNP, yNP, 'o', color = 'blue')
    plt.title('NP')
    plt.xlabel('H1N1')
    plt.ylabel('H3N2')

    plt.subplot(4,2,6)
    plt.plot(xNAb, yNAa, 'o', color = 'blue')
    plt.title('NAa/b')
    plt.xlabel('H1N1')
    plt.ylabel('H3N2')

    plt.subplot(4,2,7)
    plt.plot(xM, yM, 'o', color = 'blue')
    plt.title('M')
    plt.xlabel('H1N1')
    plt.ylabel('H3N2')

    plt.subplot(4,2,8)
    plt.plot(xNS, yNS, 'o', color = 'blue')
    plt.title('NS')
    plt.xlabel('H1N1')
    plt.ylabel('H3N2')

    plt.show()

def traverse_mgc():
    mgc = write_mgc()
    global errorInitial
    global errorFinal

    for i in range(len(mgc)):
        h1n1 = 0
        h3n2 = 0
        for j in range(len(mgc[i])):
            if j%2==0:
                h3n2 += mgc[i][j]
            else:
                h1n1 += mgc[i][j]
        if(h1n1 > h3n2):
            errorFinal += h3n2
        elif (h1n1 < h3n2):
            errorFinal += h1n1
        else:
            continue

def addlabels(x, y):
    for i in range(len(x)):
        plt.text(i, y[i], y[i], ha = 'center')

def plotHistogram():    
    data = {'Initial Error':errorInitial, 'Final Error':errorFinal}
    labels = list(data.keys())
    values = list(data.values())

    plt.bar(labels, values, color ='green', width = 0.4)
    addlabels(labels, values)
    plt.ylabel("Number of incorrect UMIs")
    plt.show()

# write_mgc()
# SumPlotGraph()
# individualInitialPlotGraph()
# individualFinalPlotGraph()
# traverse_mgc()
# plotHistogram()