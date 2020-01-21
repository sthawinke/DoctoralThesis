#' A function to estimate the dispersion parameter of the negative binomial distribution using full maximum likelihood, given the mean component
#'
#'  @param Y a vector of counts
#'  @param X a matrix of regressors
#'  @param Z a matrix of offsets
#'  @param coefs a vector of coefficients
#'  @param odOld the previous estimate of the overdispersion
#'  @param ... more arguments, passed on to the nleqslv() function
#'
#'  @return The new value of the overdispersion
#'
#' @importFrom nleqslv nleqslv
estOverdispersion = function(Y, mu, odOld, ...){
  nleqslv(x = odOld, fn = odScore, jac = NULL, mu = mu, y = Y,...)$x
}