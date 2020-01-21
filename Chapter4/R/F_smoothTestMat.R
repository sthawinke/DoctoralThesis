smoothTestMat = function(mat, x, polynoms = list(b3 = function(y,...){1+y+y^2+y^3},
                                                 b4 = function(y,...){1+y+y^2+y^3+y^4}),
                         orthoFunsList = c(h0,beta0ScoreFun,odScore), nCores = 1,
                         truncQuantLow = 1e-6, truncQuant = 1-truncQuantLow,
                         sEst = NULL, maxit = 25, prevCutOff = 1, single = TRUE,
                         xFac = NULL, maxNum = 2^20, returnEvals = FALSE,...){
if(all(mat==0)){return(matrix(NA, 3, ncol(mat),
                              dimnames = list(c(names(polynoms) ,"testStat"),
                                              colnames(mat))))}
rowID = rowSums(mat)>0
mat = mat[rowID, colMeans(mat==0)<prevCutOff, drop = FALSE]
s = if(is.null(sEst)) rowSums(mat) else sEst
x = x[rowID,]
p = NCOL(mat)
d = NCOL(x)
n = NROW(mat)
h0 = function(y, ...){rep.int(1L, length(y))}
orthLength = length(orthoFunsList)
nDisp = if(single) 1 else d
degree = length(polynoms)
polyNames = if(is.null(names(polynoms))) seq_along(polynoms) else names(polynoms)
polynoms = lapply(seq_along(polynoms), function(x){
  list(func = polynoms[[x]], norm = 1L, orthFactors = integer(orthLength+x-1))
})
vcovNB = matrix(0,d+nDisp,d+nDisp)
scoreStats = mclapply(seq_len(p), mc.cores = nCores, function(j){
  y = mat[,j]
  NBfit2 = glm.nb2(y =y, reg = x, s = s, single = single, xFac = xFac,...)
  if(anyNA(NBfit2$vcov)) {return(c(score = rep(NA, degree), testStat = NA))}
  muEst = exp(x %*% NBfit2$betas)*s
  odEst = NBfit2$theta
  vcovNB[seq_len(d),seq_len(d)] = NBfit2$vcov
  diag(vcovNB)[d+seq_len(nDisp)] = attr(NBfit2$theta,"SE")^2
  if(single) {odEst = rep(odEst, n)}

hlist = lapply(seq_along(y), function(i){
  Max = max(min(max(qnbinom(truncQuant, mu = muEst[i], size = odEst[i]),1),
                maxNum), y[i])
  ySeq = c(0,seq_len(Max))
  densVec = dnbinom(ySeq, mu = muEst[i], size = odEst[i])
  #Pre evaluate and calculate norms for the score function
  odScoreFunEval = odScore(mu = muEst[i], od = odEst[i], x = x[i,],
                           orthoFunsList = orthoFunsList, y = ySeq,
                           ySeqCum = c(0,cumsum(1/(ySeq+odEst[i]))))
  odScoreFunDensVecNorm = odScoreFunEval*densVec/sum(densVec*odScoreFunEval^2)
  #Get the parameters while still in the polynomial regime,
  #so we can use the moment generating functions
  mgf = c(1,mgfNB(mu = muEst[i], theta = odEst[i]))
  AB = odEst[i]/(odEst[i]+muEst[i])*c(-muEst[i],1)
  # The h1 polynomial coefficients
  fifthMom =  getFifthMom(mu = muEst[i], od = odEst[i])

  polyParamsVec3 = c(1 - sum(mgf[-5]),rep.int(1L,3))#h0
  #sum(tcrossprod(c(polyParamsVec3,0), y)*densVec)
  polyParamsVec3[1:2] = polyParamsVec3[1:2] - c((polyParamsVec3*AB[2]) %*%
                                                  mgf[-1])/#The first component has just been made zero
  sum(mgf[1:3]*c(AB[1]^2, 2*prod(AB), AB[2]^2))*AB #beta0ScoreFun
  #sum(tcrossprod(c(polyParamsVec3,0), yMat)*densVec*ySeq)
  polyParamsVec4 = c(1 - sum(mgf), rep.int(1L,4))#h0
  polyParamsVec4[1:2] = polyParamsVec4[1:2] -
    c((polyParamsVec4*AB[2]) %*% c(mgf[-1], fifthMom))/
    sum(mgf[1:3]*c(AB[1]^2, 2*prod(AB), AB[2]^2))*AB

  targetFunEvals = tcrossprod(cbind(1, ySeq, ySeq^2, ySeq^3, ySeq^4),
                              rbind(c(polyParamsVec3,0), polyParamsVec4))
  targetFunEvals = targetFunEvals - odScoreFunEval %o%
    c(odScoreFunDensVecNorm %*% targetFunEvals) #overdispersion score function

#Now only orthogonalize h3 and h4
  targetFunEvals[,2] = targetFunEvals[,2] -
    sum(densVec * targetFunEvals[,1] * targetFunEvals[,2])/
    sum(targetFunEvals[,1]^2*densVec)*targetFunEvals[,1]

  #Normalize
  targetFunEvals = targetFunEvals %*% (diag(sqrt(1/c(densVec %*% targetFunEvals^2))))

  #Variance
    # IetaThetaEntries = (digamma(ySeq + odEst[i])*densVec) %*% targetFunEvals
    if(returnEvals){
list(targetFunEvals = targetFunEvals, orthFunEvals = rbind(AB[1]+AB[2]*ySeq, odScoreFunEval, 1), densVec = densVec)
    } else {
      return(list(targetFunEvals = targetFunEvals))
    }
})
testStatVar = diag(degree)
score = sapply(seq_along(y), function(i){
  hlist[[i]]$targetFunEvals[y[i]+1,]
})
scoreSum = if(degree==1) sum(score) else rowSums(score)
if(returnEvals){
  list(scoreMat = c(score = scoreSum/sqrt(n),
                    testStat = scoreSum %*% testStatVar %*% scoreSum/n),
       funEvals = lapply(hlist, function(dylan){
         dylan[c("targetFunEvals", "orthFunEvals", "densVec")]}))
} else {c(score = scoreSum/sqrt(n),
          testStat = scoreSum %*% testStatVar %*% scoreSum/n)}
})
if(returnEvals){
 list(scoreStats = sapply(scoreStats, function(x){x$scoreMat}),
      funEvals = lapply(scoreStats, function(x){x$funEvals}))
} else {
  scoreStats = matrix(unlist(scoreStats), nrow = degree+1)
rownames(scoreStats)= c(polyNames ,"testStat")
colnames(scoreStats) = colnames(mat)
scoreStats
}
}
