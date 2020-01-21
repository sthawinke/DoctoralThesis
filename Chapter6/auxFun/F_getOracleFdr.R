#' Return the oracle fdr estimates, either by estimating the null density
#' directly on the null z-values, or by using weights 0 and 1 in the weighted
#' likelihood estimation
#' @param obj A fitted pransi object
#' @param nullID an indicator vector for true null values
#' @param estP0args additional arguments for the estP0 function
#'
#' @return A list with components
#' \item{f0}{An estimate of the null density}
#' \item{fdr}{The associated fdr}
#' \item{Fdr}{The associated Fdr}
getOracleFdr  = function(obj, nullID,
                         estP0args = list(z0quantRange = seq(0.05, 0.45, 0.0125),
                                          smooth.df = 3),
                         ...){
    fit = estNormal(obj$statObs[nullID]) #Fit the normal density
    p0 = do.call(estP0,
                 c(list(statObs = obj$statObs,
                        nullDensCum = pnorm(obj$zSeq, fit["mean"], fit["sd"]),
                        zSeq = obj$zSeq),
                   estP0args)) #Estimate p0
    fdrs = getFdr(obj$statObs, G0 = pnorm(obj$zSeq, fit["mean"], fit["sd"]),
                  g0 = dnorm(obj$zSeq, fit["mean"], fit["sd"]), p = length(obj$statObs),
                  p0  = p0, zValsDensObs = obj$zValsDensObs, smoothObs = TRUE,
                  zSeq = obj$zSeq, fdr = NULL)
    return(fdrs)
}
