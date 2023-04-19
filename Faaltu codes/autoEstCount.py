'''
The function autoEstCont is used to estimate the contamination fraction in a SoupChannel object, where the SoupChannel contains multiple distinct 
cell types with different marker genes. The function calculates the contamination fraction at the cluster level for each of the marker genes, and 
clusters are pruned to remove those that give implausible estimates. The function first identifies the marker genes using the tfidf method, filters
the genes by their tf-idf value and expression in the soup, and then calculates the posterior distribution of the contamination fraction for each 
cluster/gene pair. The most probable value of the contamination fraction is then taken as the final global contamination fraction. The function 
takes several parameters, such as topMarkers, tfidfMin, soupQuantile, maxMarkers, contaminationRange, rhoMaxFDR, priorRho, priorRhoStdDev, doPlot, 
forceAccept, and verbose. The function returns a modified SoupChannel object where the global contamination rate has been set and information about 
the estimation is stored in the slot fit.

The function then calculates the posterior distribution of the contamination fraction for each gene/cluster pair using a gamma prior, controlled 
by the parameters priorRho and priorRhoStdDev. These posterior distributions are then aggregated to produce a final estimate of the contamination 
fraction. Finally, the function returns a modified SoupChannel object where the global contamination rate has been set and information about the 
estimation is stored in the slot fit.

The autoEstCont() function assumes that the input SoupChannel object contains multiple distinct cell types with different marker genes. If you try 
to run it on a channel with very homogeneous cells, such as a cell line or flow-sorted cells, you may get a warning, an error, and/or an extremely 
high contamination estimate. In such cases, the function suggests manually setting the contamination to something reasonable.

The function includes several parameters to control the behavior of the algorithm, such as topMarkers, tfidfMin, soupQuantile, maxMarkers, 
contaminationRange, rhoMaxFDR, priorRho, priorRhoStdDev, doPlot, forceAccept, and verbose. These parameters allow users to adjust the specificity 
of marker genes, the minimum tf-idf value, the expression quantile in the soup, the maximum number of markers to keep, the range of acceptable 
contamination fractions, the false discovery rate, the gamma prior on the contamination fraction, and other aspects of the algorithm.
'''

def autoEstCont(sc, topMarkers=None, tfidfMin=1.0, soupQuantile=0.90, maxMarkers=100, contaminationRange=[0.01,0.8], rhoMaxFDR=0.2, priorRho=0.05, priorRhoStdDev=0.10, doPlot=True, forceAccept=False, verbose=True):
    import numpy as np
    import pandas as pd
    import estimateNonExpressingCells
    from scipy.stats import poisson
    from statsmodels.stats.multitest import multipletests
    
    if 'clusters' not in sc.metaData.columns:
        raise ValueError("Clustering information must be supplied, run setClusters first.")

    #First collapse by cluster
    s = sc.metaData.groupby('clusters')['index'].apply(list).to_dict()
    tmp = pd.DataFrame([sc.toc.loc[:,v].sum(axis=1) for k,v in s.items()], index=s.keys()).T
    ssc = sc.copy()
    ssc.toc = tmp
    ssc.metaData = pd.DataFrame({'nUMIs': tmp.sum(axis=0)}, index=tmp.columns)

    ###################
    # Get best markers
    #Get the top N soup Genes
    soupProf = ssc.soupProfile.sort_values('est', ascending=False)
    soupMin = soupProf['est'].quantile(soupQuantile)
    #Find or load markers.
    if topMarkers is None:
        #Refine this to the best markers we can manage
        mrks = quickMarkers(sc.toc, sc.metaData.clusters, N=np.inf)
        #And only the most specific entry for each gene
        mrks = mrks.sort_values(['gene', 'tfidf'], ascending=[True,False])
        mrks = mrks.drop_duplicates(subset='gene', keep='first')
        #Order by tfidif maxness
        mrks = mrks.sort_values('tfidf', ascending=False)
        #Apply tf-idf cut-off
        mrks = mrks[mrks['tfidf'] > tfidfMin]
    else:
        mrks = topMarkers
    #Filter to include only those that exist in soup 
    tgts = soupProf.loc[soupProf['est'] > soupMin].index.intersection(mrks['gene'])
    #And get the ones that pass our tfidf cut-off
    filtPass = mrks.loc[mrks['gene'].isin(tgts)].sort_values('tfidf', ascending=False)
    tgts = filtPass['gene'].head(maxMarkers).tolist()
    if verbose:
        print(f"{len(mrks)} genes passed tf-idf cut-off and {len(filtPass)} soup quantile filter.  Taking the top {len(tgts)}.")
    if len(tgts) == 0:
        raise ValueError("No plausible marker genes found.  Is the channel low complexity?  If not, reduce tfidfMin or soupQuantile")
    if len(tgts) < 10:
        warnings.warn("Fewer than 10 marker genes found.  Is this channel low complexity?  If not, consider reducing tfidfMin or soupQuantile")
        
    ############################
    # Get estimates in clusters
    #Get which ones we'd use and where with canonical method
    tmp = {k: v for k,v in zip(tgts, [1]*len(tgts))}
    ute = estimateNonExpressingCells(sc, tmp, maximumContamination=max(contaminationRange), FDR=rhoMaxFDR)
    ute = pd.DataFrame(ute, columns=sc.metaData['clusters'], index=tgts).T
    #Now calculate the observed and expected counts
    expCnts = np.outer(ssc['soupProfile']['est'], ssc['metaData']['nUMIs'])
    expCnts = pd.DataFrame(expCnts)
    expCnts.index = ssc['soupProfile'].index
    expCnts.columns = ssc['metaData'].index
    expCnts = expCnts.loc[tgts, :]

    obsCnts = ssc['toc'].loc[tgts, :]
    pp = poisson.cdf(obsCnts, expCnts*max(contaminationRange))
    qq = multipletests(pp, method='fdr_bh')[1]
    qq = pd.DataFrame(qq, index=pp.index, columns=pp.columns)

    rhos = obsCnts / expCnts
    rhoIdx = np.apply_along_axis(lambda e: np.argsort(np.argsort(e)), 1, rhos)
    rhoIdx = pd.DataFrame(rhoIdx.T, index=rhos.index, columns=rhos.columns)

    dd = pd.DataFrame({
        'gene': np.repeat(ssc['ute'].index, ssc['ute'].shape[1]).tolist(),
        'passNonExp': ssc['ute'].values.flatten(),
        'rhoEst': rhos.values.flatten(),
        'rhoIdx': rhoIdx.values.flatten(),
        'obsCnt': obsCnts.values.flatten(),
        'expCnt': expCnts.values.flatten(),
        'isExpressedFDR': qq.values.flatten(),
        'geneIdx': ssc['mrks'].reset_index().set_index('gene').loc[:, 'index'].loc[ssc['ute'].index].tolist(),
        'tfidf': ssc['mrks'].loc[ssc['ute'].index, 'tfidf'].tolist(),
        'soupIdx': ssc['soupProf'].reset_index().set_index('index').loc[ssc['ute'].index, :].index.tolist(),
        'soupExp': ssc['soupProf'].reset_index().set_index('index').loc[ssc['ute'].index, 'est'].tolist(),
    })
    dd['useEst'] = dd['passNonExp']
    if sum(dd['useEst']) < 10:
        print("Fewer than 10 independent estimates, rho estimation is likely to be unstable. Consider reducing tfidfMin or increasing SoupMin.")
    if verbose:
        print(f"Using {sum(dd['useEst'])} independent estimates of rho.")

# Last 70 lines of code left to be converted to python
'''
#Now aggregate the posterior probabilities for the ones we're including
  p.L = function(x,alpha){if(x==0){0}else{qgamma(alpha,x)}}
  p.U = function(x,alpha){qgamma(1-alpha,x+1)}
  alpha=0.95
  alpha=(1-alpha)/2
  dd$rhoHigh=sapply(seq(nrow(dd)),function(e) p.U(dd$obsCnt[e],alpha)/dd$expCnt[e])
  dd$rhoLow=sapply(seq(nrow(dd)),function(e) p.L(dd$obsCnt[e],alpha)/dd$expCnt[e])
  rhoProbes=seq(0,1,.001)
  #Using 95% confidence intervals
  #tmp = sapply(rhoProbes,function(e) {w=which(dd$useEst & dd$tfidf<1.5);sum(e>=dd$rhoLow[w] & e<=dd$rhoHigh[w])/length(w)})
  #Do a posterior estimation instead.  Use gamma prior defined by mode (priorRho) and standard deviation (priorRhoStdDev), which yields a posterior distribution for gamma of the form dgamma(rho,obsCnt+k,scale=theta/(1+theta*expCnts)). Where k and theta are the parameters for prior distribution derived using the above constraints.
  v2 = (priorRhoStdDev/priorRho)**2
  k = 1 +v2**-2/2*(1+sqrt(1+4*v2))
  theta = priorRho/(k-1)
  tmp = sapply(rhoProbes,function(e) {
                 tmp = dd[dd$useEst,]
                 mean(dgamma(e,k+tmp$obsCnt,scale=theta/(1+theta*tmp$expCnt)))
                  })
  #Calculate prior curve
  xx=dgamma(rhoProbes,k,scale=theta)
  #Get estimates
  w = which(rhoProbes>=contaminationRange[1] & rhoProbes<=contaminationRange[2])
  rhoEst = (rhoProbes[w])[which.max(tmp[w])]
  rhoFWHM = range((rhoProbes[w])[which(tmp[w]>=(max(tmp[w])/2))])
  contEst = rhoEst
  if(verbose)
    message(sprintf("Estimated global rho of %.2f",rhoEst))
  ##I think the best way to do this is based on the density.
  #tmp = density(dd$rhoEst[dd$useEst],...)
  #contEst = tmp$x[which.max(tmp$y)]
  if(doPlot){
    plot(rhoProbes,tmp,'l',
         xlim=c(0,1),
         ylim=c(0,max(c(xx,tmp))),
         frame.plot=FALSE,
         xlab='Contamination Fraction',
         ylab='Probability Density')
    #Add prior
    lines(rhoProbes,xx,lty=2)
    abline(v=rhoProbes[which.max(tmp)],col='red')
    legend(x='topright',
           legend=c(sprintf('prior rho %g(+/-%g)',priorRho,priorRhoStdDev),
                    sprintf('post rho %g(%g,%g)',rhoEst,rhoFWHM[1],rhoFWHM[2]),
                    'rho max'),
           lty=c(2,1,1),
           col=c('black','black','red'),
           bty='n')
    #plot(0,
    #     xlim=c(0,1),
    #     ylim=c(0,max(tmp$y)),
    #     type='n',
    #     frame.plot=FALSE,
    #     xlab='Contamination Fraction',
    #     ylab='Density'
    #     )
    #lines(tmp$x,tmp$y)
    #abline(v=contEst,col='red')
  }
  sc$fit = list(dd=dd,
                priorRho=priorRho,
                priorRhoStdDev=priorRhoStdDev,
                posterior = tmp,
                rhoEst = rhoEst,
                rhoFWHM = rhoFWHM,
                markersUsed = mrks
                )
  #Set the contamination fraction
  sc = setContaminationFraction(sc,contEst,forceAccept=forceAccept)
  return(sc)
}
'''