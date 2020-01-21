# A function to calculate the ratio between the inner product of a target and orthogonalizing function and the norm of the orthogonalizing function, all in a numeric way
getInnerScoreRatio = function(targetFunEval, orthFunDensVecNorm){
  orthFunDensVecNorm %*% targetFunEval
}
