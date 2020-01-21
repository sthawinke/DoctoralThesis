#' A function to generate NB data
genNB = function(NBcoef, NBtheta, x, libSizes, idDA, beta,
                 mu = rowMultiply(exp(NBcoef + outer(idDA, x)*beta), libSizes),
                 Sigma = NULL,   n = length(libSizes), p = length(NBcoef),
                 taxNames = seq_len(p), samNames = names(libSizes),
                 empirical = FALSE){
mat = if(is.null(Sigma)) {
  matrix(rnbinom(n*p,size = NBtheta, mu = mu), n, p, byrow = TRUE)
} else {
    Rmvnegbin(n = n, mu = mu, Sigma = Sigma, ks = NBtheta, empirical = empirical)
}
taxNames[as.logical(idDA)] = paste0(taxNames[as.logical(idDA)], "-TP")
dimnames(mat) = list(samNames, taxNames)
mat
}
genNBsim = function(Cor, p, ...){
  args = list(...)
  if(Cor %in% c("None","Cor")){
    args$FC = 1
  }
taxaID = order(args$coef, decreasing = TRUE)[seq_len(p)]
args$coef = args$coef[taxaID]
args$coef = args$coef - log(sum(exp(args$coef)))
args$theta = args$theta[taxaID]
  do.call(genNB, args)
}