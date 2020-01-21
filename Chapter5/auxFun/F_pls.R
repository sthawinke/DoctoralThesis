#PLS
plsSim = function(datList, what = rep("pca", length(datList)), dim = 3L){
    matList = transformData(datList, what)
    fit = block.pls(X = matList, indY = 1L, ncomp = dim, mode = "canonical")
    for (i in seq_along(datList)){
        colnames(fit$variates[[i]]) = colnames(fit$loadings[[i]]) =
            paste0("Dim", seq_len(dim))
    }
    list("featureLoadings" = fit$loadings, "sampleScores" = fit$variates)
}