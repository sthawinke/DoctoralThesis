#' A function to obtain a weighted inner product of two functions, weighted by the negative binomial density.
getInnerProduct = function(fun1, fun2, mu, od, x, ySeq = ySeq, ySeqCum = ySeqCum,densVec,...){
  sum((if(class(fun1)=="function") fun1(y = ySeq, mu = mu, od = od, x = x, ySeqCum = ySeqCum,...) else evalTargetFun(fun1, y = ySeq,mu = mu, od = od, x = x, ySeqCum = ySeqCum,...)) *(if(class(fun2)=="function") fun2(y = ySeq, mu = mu, od = od, x = x, ySeqCum = ySeqCum,...) else evalTargetFun(fun2, y = ySeq, mu = mu, od = od, x = x, ySeqCum = ySeqCum,...))* densVec)
}