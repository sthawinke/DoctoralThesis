#' A function to calculate the variance of the test statistic.
#'
testStatisticVar = function(mu, od, x, NBvcov, polynoms, degree = length(polynoms[[1]]), diGamma, truncQuant = 1-1e-8){
  IetaTheta = matrix(rep(0, 3*degree), ncol = degree)
  IetaTheta[nrow(IetaTheta),] = vapply(seq_along(polynoms[[1]]), FUN.VALUE = 0, function(z){
    sum(vapply(seq_along(mu), FUN.VALUE = 0, function(i){getCovariance(match.fun(polynoms[[i]][[z]]), diGamma, mu = mu[i], od = od, x = x[i], truncQuant = truncQuant)}))
  })
solve((diag(rep(length(mu), degree), nrow =degree, ncol = degree) - crossprod(IetaTheta, NBvcov)  %*%  IetaTheta)/length(mu))
}
