#' A function that evaluates the Jacobian for the estimation of the overdispersion parameter
#'
#' @param x the current value of the estimate
#' @param Y a vector of counts
#' @param mu a vector of means
#' @param truncQuant the quantile at which to truncate
#'
#' @return a scalar, the evaluation of the jacobian function at x
odJac = function(x, y, mu, truncQuant = 0.99, ...){
  tmp = sapply(mu, function(mui){
    truncAtSeq = c(0, seq_len(qnbinom(p = truncQuant, mu = mui, size = 1/x)))
    sum(1/(1/x+ truncAtSeq)^2 * (1-pnbinom(size = 1/x, mu = mui, q = truncAtSeq, lower.tail = TRUE)))})
  sum(1/x^4*(tmp - x*mu/(mu+1/x)))
  }