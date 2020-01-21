# JIVE
# See Lock et al., 2013, OConnell and Lock 2018
jiveSim = function(datList, method = "given", dim = 3, jointDim = 0L, ...){
    n_indiv = rep(dim, length(datList))
    result = jive(lapply(extractData(datList), t), method = method,
                  rankJ = jointDim, rankA = n_indiv, showProgress = FALSE,...)
    l <- length(result$data)
    nPCs = sum(n_indiv)
    PCs = list("matrix", length(datList))
    for (i in seq_len(l)) {
        if (n_indiv[i] > 0) {
            SVD = svd(result$individual[[i]], nu = n_indiv[i],
                      nv = n_indiv[i])
            PCs[[i]] = t(diag(SVD$d)[seq_len(n_indiv[i]), seq_len(n_indiv[i])] %*%
                    t(SVD$v[, seq_len(n_indiv[i])]))
            colnames(PCs[[i]]) = paste0("Dim", seq_len(n_indiv[i]))
            rownames(PCs[[i]]) = rownames(datList[[i]])
            }
        }
    names(PCs) = names(datList)
    list("sampleScores" = PCs, "converged" = result$converged)
}