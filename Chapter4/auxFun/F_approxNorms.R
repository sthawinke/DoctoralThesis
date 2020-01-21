#' A function to approximate the norms of the score functions numerically, for a single observation and corresponding parameters
#'
#' @param beta the vector of regression coefficients for the mean
#' @param X the vector of regressors
#' @param Z the vector of offsets
#' @param od the overdispersion parameter
#' @param truncQuant the quantile at which to truncate the approximation
#'
#' @return a vector with all the score norms
approxNorms = function(beta, X, Z, od, truncQuant = 1-1e-15){
mu = exp(sum(X * beta) + sum(Z))
ySeq = c(0L, seq_len(qnbinom(p = truncQuant, mu = mu, size = od)))
densVec = dnbinom(ySeq, mu = mu, size =od)
scoreNormBeta0 = sum(((ySeq-mu)/(1+mu/od))^2*densVec)
scoreNormBeta1 = sum((X[-1]*(ySeq-mu)/(1+mu/od))^2*densVec)
ySeq2 = ySeq[-c(1:2)]-1
tmp = c(0,0,cumsum(ySeq2/(1+ySeq2/od)))
scoreNormod = sum((od^2*log(1+mu/od) -(ySeq+od)*mu/(1+mu/od) + tmp)^2*densVec)
c(beta0 = scoreNormBeta0, beta1 = scoreNormBeta1, od = scoreNormod)
}