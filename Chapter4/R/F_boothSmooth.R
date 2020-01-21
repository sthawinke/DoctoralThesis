#' A bootstrapping function for the smooth test statistic.
#'
#' @details The library sizes are assumed known and passed on the the smooth testing function, rather than be recalulated
bootSmooth = function(mu, ods, B, Libs,x,...){
  n = length(mu)
  dataMat = matrix(rnbinom(n = n*B, size = ods, mu = mu), n,B, byrow = FALSE)
  smoothTestMat(dataMat, x, sEst = Libs,...)["testStat",]
}