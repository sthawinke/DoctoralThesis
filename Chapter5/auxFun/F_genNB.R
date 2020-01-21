#' Generate negative binomial data
genNB = function(estList, nPop, n, p, FC, TPR, compensation = TRUE){
    if((n %% nPop)!=0) stop("Choose even group sizes!")
    groupIndComp = rep(seq_len(nPop), times = n/nPop)
    p = min(p, length(estList$coef))
    ## Sample parameters
    rhosSampled0 = exp(sample(estList$coef, p))
    rhosSampled0 = rhosSampled0/sum(rhosSampled0)
    #   For each population, modify a fraction TPR of the abundances (compensation is not so badly needed here, we do not test     individual species). Record the taxa that have been modified.
    rhosIDSampled = lapply(integer(nPop), function(x){
        addFoldChange(rhosSampled0, fc = FC, H1frac = TPR, compensate = compensation)
    })
    rhosSampled = simplify2array(rhosIDSampled)[,groupIndComp]
    thetasSampled = sample(estList$theta, p)
    libSizesSampled = sample(estList$Libs, n) #Sampled libSizes
    ## Create dataset
    meanMat = libSizesSampled * t(rhosSampled)
    dataMat = matrix(rnbinom(n*p, mu = meanMat, size = matrix(thetasSampled,n,p,byrow = TRUE)),n,p)
    colnames(dataMat) = rownames(rhosSampled)
    rownames(dataMat) = paste0("Group", groupIndComp, "_", seq_len(n))
    rowID = rowSums(dataMat)>0
    colID = colSums(dataMat)>0
    dataMat = dataMat[rowID, colID]
    list(dataMat = dataMat, libSizesSampled = libSizesSampled[rowID],
         thetasSampled = thetasSampled[colID], rhosSampled = lapply(rhosIDSampled, function(x) x[colID]))
}