#' The score function for the regression parameters of the negative binomial distribution
#'
#' @param y the observed count
#' @param mu the expected count
#' @param od the overdispersion
#' @param x the regressor

betaScoreFun = function(y, mu, od, x,...){x*(y-mu)/(1+mu/od)}
beta0ScoreFun = function(y, mu, od,x,...){(y-mu)/(1+mu/od)}

