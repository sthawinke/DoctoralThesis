#' A function to approximate the inner products of the score functions numerically
#'
#' @param beta the vector of regression coefficients for the mean
#' @param X the vector of regressors
#' @param Z the vector of offsets
#' @param od the overdispersion parameter
#' @param truncQuant the quantile at which to truncate the approximation
#' @param coef vector of coefficients of the polynomial
#'
#' @return a vector with all the inner products
approxInnerProducts = function(beta, X, Z, od, coefs, truncQuant = 1-1e-15){
  mu = exp(sum(X * beta) + sum(Z))
  ySeq = c(0L, seq_len(qnbinom(p = truncQuant, mu = mu, size = od)))
  densVec = dnbinom(ySeq, mu = mu, size =od)
  degree = length(coefs)-1
  polyN = t(coefs*t(model.matrix(data = data.frame(ySeq = ySeq),object = formula(paste0("~",paste0("I(ySeq^",1:degree,")", collapse ="+"))))))
  innnerProdBeta0 = sum(rowSums(((ySeq-mu)/(1+mu/od))*polyN)*densVec)
  innnerProdBeta1 = sum(rowSums((X[2]*(ySeq-mu)/(1+mu/od))*polyN)*densVec)
  ySeq2 = ySeq[-c(1:2)]-1
  tmp = c(0,0,cumsum(ySeq2/(1+ySeq2/od)))
  innnerProdod = sum(rowSums((od^2*log(1+mu/od) -(ySeq+od)*mu/(1+mu/od) + tmp)*polyN)*densVec)
  c(beta0 = innnerProdBeta0, beta1 = innnerProdBeta1, od = innnerProdod)
}