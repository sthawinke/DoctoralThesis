#' Combine all other methods in one function
comboSim = function(datList, dim = 3, clrT = rep(TRUE, length(datList))){
    l = length(datList)
    #Prune
    datList = lapply(datList, function(x){x[, colSums(x)>0]})
    ## Multiple CCA, no shrinkage
    ccaFitUnshrunk = ccaSim(datList)
    ## Multiple CCA, with shrinkage
    ccaFitShrunk = ccaSim(datList, penalty = TRUE)
    ## Multiple CCA, on clr tranformed data
    ccaFitUnshrunkClr = ccaSim(datList, clrT = clrT)
    ## Jive
    jiveFit = jiveSim(datList)
    ## MoFa
    mofaFit = mofaSim(datList, maxIter = 1e3L)
    ## PCA, scaled
    pcaFit = pcaSim(datList, what = rep("pca", l))
    ## Correspondence analysis
    correspFit = pcaSim(datList, what = rep("ca", l))
    ## PCA, clr
    pcaFitClr = pcaSim(datList, what = ifelse(clrT, "clr", "ca"))
    ## PLS
    plsFit = plsSim(datList)
    ## PLS clr
    plsFitClr = plsSim(datList, what = ifelse(clrT, "clr", "pca"))

    list("ccaUnshrunk" = ccaFitUnshrunk, "ccaShrunk" = ccaFitShrunk,
         "ccaUnshrunkClr" = ccaFitUnshrunkClr, "jive" = jiveFit,# "mofa" = mofaFit,
         "pca" = pcaFit, "correspFit" = correspFit, "pcaFitClr" = pcaFitClr,
         "pls" = plsFit, "plsClr" = plsFitClr, "mofaFit" = mofaFit)
}
