#' Functions to generate multivariate data and a covariance matrix.
makeZnorm = function(n, p, p0, FC, Sigma, empirical = FALSE,
                     oneGroup = FALSE){
  mu2 = c(rep(0, p*p0), rep(FC, p*(1-p0)))

  if(oneGroup){
  normMat = SpiecEasi::rmvnorm(n = n, Sigma = Sigma, mu = mu2,
                                  empirical = empirical)
  } else {
  normMat1 = SpiecEasi::rmvnorm(n = n/2, Sigma = Sigma, mu = integer(p),
                                empirical = empirical)
  normMat2 = SpiecEasi::rmvnorm(n = n/2, Sigma = Sigma, mu = mu2,
                                empirical = empirical)
  normMat = rbind(normMat1, normMat2)
  }
  colnames(normMat) = seq_len(p)
  colnames(normMat)[mu2!=0] = paste(colnames(normMat)[mu2!=0], "-TP")
  return(normMat)
}
makeSigma = function(p, Sd = 1, tol = 1e-6, random = FALSE, lower = -1, upper = 1){
  if(!random){
    tmp = matrix(Sd, p, p)
    diag(tmp) = 1
    return(tmp)
  } else {
    ranVec = rnorm(p, sd = Sd)
    tmp = outer(ranVec, ranVec)
    # if((lower < -1) | (upper > 1)) stop("Correlations must be in [-1,1] range")
    # tmp = matrix(0, p, p)
    # tmp[upper.tri(tmp)] = tmp[lower.tri(tmp)] = runif((p^2-p)/2, min = lower, max = upper)
    diag(tmp) = 1
    ev = eigen(tmp, symmetric = TRUE)$values
    if(!all(ev >= -tol * abs(ev[1L]))){
      tmp = makeSigma(p, Sd = Sd, random = random, lower = lower, upper = upper)
    }
    return(tmp)
  }
}