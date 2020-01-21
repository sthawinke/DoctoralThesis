#' Check for convergence of iterative data integrations
meanMOFAconv = function(resList){
    mean(sapply(resList, function(x) x$mofaFit$converged))
}
meanCompIntconv = function(resList){
    tmp = rowMeans(sapply(resList, function(x) x$converged))
    names(tmp) = paste0("Dim", seq_along(tmp))
    tmp
}