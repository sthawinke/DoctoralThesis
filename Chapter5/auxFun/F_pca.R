# PCA/CA
pcaSim = function(datList, dim = 2, what = rep("pca", length(datList))){
    matList = transformData(datList, what)
    concat = Reduce(matList, f = cbind)
    Svd = svd(concat)
    # Extract sample scores
    sampleScores = (Svd$u %*% diag(Svd$d))[, seq_len(dim)]
    rownames(sampleScores) = rownames(concat)
    # Extract feature loadings
    featureLoadings = Svd$v[, seq_len(dim)]
    colnames(featureLoadings) = colnames(sampleScores) = paste0("Dim", seq_len(dim))
    #Divide into a list again
    featureLoadings = lapply(datList, function(dat){
        tmp = featureLoadings[seq_len(ncol(dat)),]
        rownames(tmp) = colnames(dat)
        tmp})
    list("featureLoadings" = featureLoadings,
         "sampleScores" = sampleScores)
}