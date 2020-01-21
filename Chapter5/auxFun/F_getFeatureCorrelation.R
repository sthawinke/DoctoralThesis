#' Calculate correlations between all features and DA features
#' @param simSeqRes The model integration list
#' @param simSeqDat The data list
#' @param simSeqTestResList A list of test results
#'
#' @return The wilcoxon rank sum test statistics converted to z-values
getFeatureCorrelation = function(simSeqRes, simSeqDat, simSeqTestResList = NULL,
                                 signMat = NULL, DAfeat = NULL, stat = "wilcox"){
    SimSeqResCols = simSeqRes$cols
    gm = function(x) exp(mean(log(x)))
    diffFun = function(all, diff, stat){
        if(length(diff) > 1)
            switch(stat,
                   "wilcox" = qnorm(quantCorrect(wilcox.test(
                       diff, all, alternative = "less")$p.value)),
                   "tTest" = qnorm(quantCorrect(t.test(
                       x = diff, y = all, alternative = "less")$p.value))
                   ) else NA
        }
    #Calculate all pairwise inner products
    allPairWise = lapply(SimSeqResCols, function(li){
        collapsedMat = Reduce(li, f = rbind)
        foo = tcrossprod(collapsedMat)
        diag(foo) = 0
        foo
    })
    # All pairwise inner products between views
    allPairWiseBetween = lapply(SimSeqResCols, function(li){
        foo = apply(li[[1]], 1, "tcrossprod", li[[2]])
        rownames(foo) = rownames(li[[2]])
        foo
    })
    if(is.null(DAfeat)){
        #Get the DA taxa
        DAfeat = lapply(simSeqDat, function(x){x$DEtaxa})
    }
    allFeat = Reduce(DAfeat, f = c)
    if(is.null(signMat)){
        #Extract the pairwise correlation of differential features,
        #and correct for the signs of the differential abundance

        signMat0 = crossprod(upDownToNumeric(cbind(simSeqTestResList[[1]]$result,
                                               simSeqTestResList[[2]]$result))
                         )[allFeat,allFeat]
        signMat = signMat0 == nrow(simSeqTestResList[[1]]$result)
    }
    #Perfect correspondence over groups
    difPairWise = lapply(allPairWise, function(pairMat){
        id = rownames(pairMat) %in% allFeat
        difMat = pairMat[id, id]
        idSign = rownames(signMat) %in% colnames(difMat)
        difMatCor = difMat * t(signMat[idSign,idSign])
        difMatCor[difMatCor!=0]
    })
    resAll = mapply(allPairWise, difPairWise, FUN = diffFun, stat = stat)
    difPairWiseBetween = lapply(allPairWise, function(pairMat){
        idRow = rownames(pairMat) %in% DAfeat[[2]]
        idCol = colnames(pairMat) %in% DAfeat[[1]]
        difMat = pairMat[idRow, idCol]
        difMatCor = difMat * t(signMat[rownames(signMat) %in% colnames(difMat),
                                       colnames(signMat) %in% rownames(difMat)])
        difMatCor[difMatCor!=0]
    })
    resBetween = mapply(allPairWise, difPairWiseBetween, FUN = diffFun, stat = stat)
    list("all" = resAll, "between" = resBetween)
}