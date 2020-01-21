# GCCA and RGCCA
# Set tau to default for no shrinkage
#' @return The feature loadings
ccaSim = function(datList, dim = 3, penalty = FALSE, standardize = TRUE,
                  clrT = rep(FALSE, length(datList)), ...){
    datList = lapply(seq_along(datList), function(i){
        if(clrT[i]) clrTransformMat(datList[[i]]) else datList[[i]]
    })
    if(penalty){
        perm.out = MultiCCA.permute(xlist = datList, standardize = standardize,
                                    trace = FALSE, ...)
    }
    fit = MultiCCA(xlist = datList, penalty = if(penalty) perm.out$bestpenalties else
        sqrt(sapply(datList, ncol)), ncomponents = dim,
             standardize = standardize, trace = FALSE, ...)
    #Assign names, and create sample scores
    sampScores = vector("list", length(datList))
    for(i in seq_along(datList)){
        rownames(fit$ws[[i]]) = colnames(datList[[i]])
        colnames(fit$ws[[i]]) = paste0("Dim", seq_len(dim))
        sampScores[[i]] = scale(datList[[i]]) %*% fit$ws[[i]]
    }
    list("featureLoadings" = fit$ws, "sampleScores" = sampScores)
}