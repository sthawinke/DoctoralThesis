#' Generate normal data
genNormal = function(estList, nPop, n, p, FC, TPR){
    if((n %% nPop)!=0) stop("Choose even group sizes!")
    groupIndComp = rep(seq_len(nPop), times = n/nPop)
    p = min(p, length(estList["coef",]))
    ## Sample parameters
    coefSampled0 = sample(estList["coef",], p)
    nOTUs = round(p*TPR) #DA taxa
    meanSampled0 = lapply(integer(nPop), function(x){
        OTUids = sample(names(coefSampled0), nOTUs, replace = FALSE)
        coefSampled0[OTUids] = coefSampled0[OTUids]+FC # Add fold change up
        indTP <- names(coefSampled0) %in% OTUids
        newTaxaNames <- paste0(names(coefSampled0)[indTP], "-TPup")
        names(coefSampled0)[indTP] <- newTaxaNames
        coefSampled0
    })
    meanSampled = simplify2array(meanSampled0)[,groupIndComp]
    dataMat = t(matrix(rnorm(n*p, mean = meanSampled, sd = estList["sd",]), ncol = n, nrow = p))
    colnames(dataMat) = rownames(meanSampled)
    rownames(dataMat) = paste0("Group", groupIndComp, "_", seq_len(n))
    list(dataMat = dataMat, meanSampled = meanSampled0)
}