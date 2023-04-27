'''
This project has been completed by the following students of BITS Pilani, Hyderabad Campus:
1. Suvigya Sharma (2020A7PS0140H)
2. Mohit Agarwal (2020A7PS0189H)
'''

import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt

# taking input from csv file
tod = pd.read_csv("data.csv", index_col=0)

# seperating the empty droplets from the dataset based on a cutoff value of 20
empty_drops = tod.loc[tod['sum']<20]

# seperating the non-empty droplets from the dataset based on a cutoff value of 20
toc = tod.loc[tod['sum']>20]

# keeping the barcode column in a seperate variable for attachment in the final step
barcode = toc['Barcode']

# dropping the barcode column from the dataset
toc.drop('Barcode', axis=1, inplace=True)

# initialization of the global variables to be used in plotting of graphs 
xInitialSum = []
yInitialSum = []
errorInitial = 0
errorFinal = 0

# auxiliary function to initialize the variables to be used in the further program
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

# function to create the background expression of all genes using the empty droplets
def create_bg_gene_empty():
    geneColumns = auxiliary()[0]
    #sum of an individual gene in all empty drops , size = number of genes
    geneSum_empty = np.sum(empty_drops[geneColumns], axis=0)    # geneSum_empty is list
    #sum of all data in empty drops
    totalSum_empty = np.sum(geneSum_empty)                      # totalSum_empty is a number
    #fraction of background expression of all genes ,size = number of genes
    bg_gene_empty = geneSum_empty/totalSum_empty                # bg_gene_empty (for empty droplets) is list
    return bg_gene_empty

# function to create the background contamination fraction rho for each cell using the non-empty droplets
def create_rho_c():
    
    # initialise variables using auxiliary function
    detected = auxiliary()[1]
    h1n1_columns = auxiliary()[2]
    h3n2_columns = auxiliary()[3]
    h1n1_column_names = auxiliary()[4]
    h3n2_column_names = auxiliary()[5]
    bg_gene_empty = create_bg_gene_empty()

    # initialise global variables
    global xInitialSum, yInitialSum, errorInitial, errorFinal

    # create rho_c table to store the contamination fraction rho for each cell
    rho_c = []

    # traversing through each cell
    for i in range(len(detected)):

        # store sum of h1n1 and h3n2 UMIs in variables to compare and select primary cell
        h1n1 = detected['subsets_H1N1_sum'][i]
        h3n2 = detected['subsets_H3N2_sum'][i]

        # to be used in SumPlotGraph plotting
        xInitialSum.append(h1n1)
        yInitialSum.append(h3n2)

        # initialise local variables
        ngc_sum = 0
        nc = 0
        bgc_sum = 0

        # if h1n1 UMIs lesser than h3n2 UMIs, h3n2 is primary cell
        if(h1n1 < h3n2):
            # add h1n1 UMIs to errorInitial for using in PlotHistogram
            errorInitial += h1n1
            # ngc is UMIs that shouldn't be present in primary cell
            ngc_sum = h1n1
            # nc is total UMIs in the cell
            nc = h3n2 + h1n1

            # if h1n1 UMIs is greater than 0, then check for each gene in h1n1_columns
            if(h1n1>0):
                for k in range(len(h1n1_column_names)):
                    if(h1n1_columns[h1n1_column_names[k]][i])>0:
                        # bgc is UMIs in the cell for genes that shouldn't be present in background
                        bgc_sum += bg_gene_empty[2*k+1]
                
                # mathematically, rho_c for the cell = ngc_sum/nc*bgc_sum
                rho_c.append(ngc_sum/(nc*bgc_sum))
            else:
                # if h1n1 UMIs is 0, then rho_c for the cell = 0
                rho_c.append(0)

        # if h3n2 UMIs lesser than h1n1 UMIs, h1n1 is primary cell
        elif(h1n1 > h3n2):
            # add h3n2 UMIs to errorInitial for using in PlotHistogram
            errorInitial += h3n2
            # ngc is UMIs that shouldn't be present in primary cell
            ngc_sum = h3n2
            # nc is total UMIs in the cell
            nc = h3n2 + h1n1

            # if h3n2 UMIs is greater than 0, then check for each gene in h3n2_columns
            if(h3n2>0):
                for k in range(len(h3n2_column_names)):
                    if(h3n2_columns[h3n2_column_names[k]][i])>0:
                        # bgc is UMIs in the cell for genes that shouldn't be present in background
                        bgc_sum += bg_gene_empty[2*k]
                
                # mathematically, rho_c for the cell = ngc_sum/nc*bgc_sum
                rho_c.append(ngc_sum/(nc*bgc_sum))
            else:
                # if h3n2 UMIs is 0, then rho_c for the cell = 0
                rho_c.append(0)
        
        else:
            # if h1n1 UMIs is equal to h3n2 UMIs, then rho_c for the cell = 0
            rho_c.append(0)

    return rho_c

# function to create the observed gene count ogc for each cell using the non-empty droplets
def create_ogc():
    # variables initialised using previous functions
    rho_c = create_rho_c()
    bg_gene_empty = create_bg_gene_empty()

    # nc_total is column containing sum of UMIs in each cell
    nc_total = toc['sum']
    # specify length of rows and columns
    ogc_rows, ogc_cols = (len(rho_c), len(bg_gene_empty))
    # initialise ogc table as 0 to store the observed gene count ogc for each cell
    ogc = [[0 for i in range(ogc_cols)] for j in range(ogc_rows)]

    # traversing through each cell
    for i in range(len(rho_c)):
        for j in range(len(bg_gene_empty)):
            # mathematically, ogc for each cell = rho_c*bg_gene_empty*nc_total
            ogc[i][j] = rho_c[i]*bg_gene_empty[j]*nc_total[i]

    return ogc

# function to create the endogenous gene count mgc for each cell using the non-empty droplets
def create_mgc(toc):
    # variables initialised using previous functions
    bg_gene_empty = create_bg_gene_empty()
    ogc = create_ogc()
    geneColumns = auxiliary()[0]
    ogc_rows, ogc_cols = (len(ogc), len(ogc[0]))
    ngc = toc[geneColumns]

    # creating mgc table intialised as 0 for all cells
    mgc = [[0 for i in range(ogc_cols)] for j in range(ogc_rows)]
    # traversing through each cell
    for i in range(len(ngc)):
        for j in range(len(ngc.columns)):
            # mathematically, mgc for each cell = ngc - ogc
            mgc[i][j] = ngc.iloc[i][j] - ogc[i][j]
            # if mgc is negative, then it is set to 0
            if(mgc[i][j]<0):
                mgc[i][j] = 0
    return mgc

# function to write the mgc table to a csv file
def write_mgc():
    # variables initialised using previous functions
    old_mgc = create_mgc(toc)
    # new_mgc created to add the barcode column and headers row
    new_mgc = []
    headers = ['Barcode', 'swH3N2_PB2_comm_raw', 'swH1N1_PB2_comm_raw', 'swH3N2_PB1_comm2_raw', 'swH1N1_PB1_comm2_raw',
                'swH3N2_PA_comm_raw', 'swH1N1_PA_comm_raw', 'swH3N2_HAa_raw', 'swH1N1_HAa_raw', 
                'swH3N2_NP_comm_raw', 'swH1N1_NP_comm_raw', 'swH3N2_NAa_raw', 'swH1N1_NAb_raw',
                'swH3N2_M_comm_raw', 'swH1N1_M_comm_raw', 'swH3N2_NS_comm_raw', 'swH1N1_NS_comm_raw']
    # headers row added to new_mgc
    new_mgc.append(headers)

    # loop to add barcode column to new_mgc for row in old_mgc
    for i in range(len(old_mgc)):
        temp = []
        temp.append(barcode[i])
        for j in range(len(old_mgc[0])):
            temp.append(old_mgc[i][j])
        new_mgc.append(temp)
    
    # new_mgc converted to dataframe and written to csv file
    mgc = pd.DataFrame(new_mgc)
    mgc.to_csv('output.csv', index=False, header=False)

    return old_mgc

# function to plot the graph which shows the initial and final total number of UMIs for each cell
def SumPlotGraph():
    # initialise variables using previous functions
    mgc = write_mgc()
    geneColumns = auxiliary()[0]
    # create variables to store final total number of UMIs in mgc table
    xFinalSum = []
    yFinalSum = []

    # loop to calculate the final total number of UMIs in mgc table
    for i in range(len(mgc)):
        h1n1 = 0
        h3n2 = 0
        for j in range (len(geneColumns)):
            # j%2==0 indicates h3n2 columns and j%2==1 indicates h1n1 columns
            if j%2==0 :
                h3n2 += mgc[i][j]
            else:
                h1n1 += mgc[i][j]
        
        xFinalSum.append(h1n1)
        yFinalSum.append(h3n2)

    # plotting the subplots in one graph
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

    # showing the plot
    plt.show()

# function to plot the graph which shows the initial number of UMIs for all genes in each cell
def individualInitialPlotGraph():
    # this mgc variable is actually toc modified, we use mgc variable name in order to maintain consistency with the finalPlotGraph function
    mgc = toc[['swH3N2_PB2_comm_raw', 'swH1N1_PB2_comm_raw', 'swH3N2_PB1_comm2_raw', 'swH1N1_PB1_comm2_raw',
                'swH3N2_PA_comm_raw', 'swH1N1_PA_comm_raw', 'swH3N2_HAa_raw', 'swH1N1_HAa_raw', 
                'swH3N2_NP_comm_raw', 'swH1N1_NP_comm_raw', 'swH3N2_NAa_raw', 'swH1N1_NAb_raw',
                'swH3N2_M_comm_raw', 'swH1N1_M_comm_raw', 'swH3N2_NS_comm_raw', 'swH1N1_NS_comm_raw']]

    # x is H1N1, y is H3N2
    # initialising lists for plotting individual genes
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

    # loop to add values to the lists
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

    # plotting the subplots in one graph
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

    # showing the plot
    plt.show()

# function to plot the graph which shows the final number of UMIs for all genes in each cell
def individualFinalPlotGraph():
    mgc = write_mgc()
    # x is H1N1, y is H3N2
    # initialising lists for plotting individual genes
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

    # loop to add values to the lists
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

    # plotting the subplots in one graph
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

    # showing the plot
    plt.show()

# function to traverse mgc to calculate final error
def traverse_mgc():
    # initialising global variables
    mgc = write_mgc()
    global errorInitial
    global errorFinal

    # calculating final error by looping through all the cells
    for i in range(len(mgc)):

        # h1n1 and h3n2 are the final number of UMIs for each gene in h1n1 and h3n2 respectively for the cell i
        h1n1 = 0
        h3n2 = 0

        # finding values of h1n1 and h3n2 by adding the values of UMIs for each gene
        for j in range(len(mgc[i])):
            if j%2==0:
                h3n2 += mgc[i][j]
            else:
                h1n1 += mgc[i][j]
        
        # adding the smaller value to the final error
        if(h1n1 > h3n2):
            errorFinal += h3n2
        elif (h1n1 < h3n2):
            errorFinal += h1n1
        else:
            continue

# function to add labels to plotted histograms
def addlabels(x, y):
    for i in range(len(x)):
        plt.text(i, y[i], y[i], ha = 'center')

# function to plot the histogram
def plotHistogram():
    # making a key-value pair for the data
    data = {'Initial Error':errorInitial, 'Final Error':errorFinal}
    labels = list(data.keys())
    values = list(data.values())

    # plotting the histogram
    plt.bar(labels, values, color ='green', width = 0.4)
    addlabels(labels, values)
    plt.ylabel("Number of incorrect UMIs")

    # showing the plot
    plt.show()

# uncomment the functions "write_mgc()" to run the code and output to output.csv
# write_mgc()

# uncomment the function below whichever plot you need to see
# SumPlotGraph()
# individualInitialPlotGraph()
# individualFinalPlotGraph()

# uncomment both the functions to get the histogram plot
# traverse_mgc()
# plotHistogram()