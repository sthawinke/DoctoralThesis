# A function to obtain Fdr and fdr
getTheorFdr = function(statObs, densFun = "dnorm",
                       zValsDensObs = do.call(densFun, c(list(x = statObs), testPargs)),
                       zSeq = NULL, testPargs = list(),
                       distFun = "pnorm", z0Quant = pnorm(c(-1,1)),
                       quantileFun = "qnorm", ...){
  if(length(z0Quant)==1) {z0Quant = sort(c(z0Quant, 1-z0Quant))}
xCutOff = do.call(quantileFun, c(list(p = z0Quant), testPargs))
ecdfVals = do.call(distFun, c(list(q = statObs), testPargs))
ecdfVals[ecdfVals>0.5] = 1-ecdfVals[ecdfVals>0.5]
steps = ecdf(statObs)(statObs)
steps[steps>0.5] = 1-steps[steps>0.5]+1/p
p00Hat = min(1,mean(statObs > xCutOff[1] & statObs < xCutOff[2])/
               diff(z0Quant))
FdrTheor = ecdfVals/steps*p00Hat
FdrTheor[FdrTheor>1] = 1
if(densFun=="dwilcox") statObs = round(statObs)
fdrTheor = zValsDensObs/approx(y = zValsDensObs, x = zSeq, xout = statObs)$y *
                     p00Hat
fdrTheor[fdrTheor>1] = 1
return(list(Fdr = FdrTheor, fdr = fdrTheor))
}