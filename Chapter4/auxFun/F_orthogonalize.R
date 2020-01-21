#' Orthogonalization functions
#'
#' @param mu a mean vector
#' @param od an overdispersion
#' @param x a vector of regressors
#' @param targetFun a character string with function name or function that needs to be orthogonalized
#' @param orthoFun a character string with function name or function with respect to which the targetFun needs to be orthogonalized
#' @param targetFunsList a character vector of functions. All are orthogonalized versus the functions in orthoFunsList, and versus the preceding functions in the list
#' @param orthoFunsList a character vector of functions to be orthogonalized against
#' @param normalize a boolean, should final functions be normalized?
#'
#' @return an orthogonalized function with arguments mu, od, x and y
#'
#' @details
orthogonalizeList = function(orthFunDensVecNorm, densVec, orthFunEvals, normalize = TRUE, orthLength = length(orthFunEvals), targetFunEvals, targetFunsList,...){
for (i in seq_along(targetFunsList)){
  for (j in seq_len(orthLength)){
    targetFunsList[[i]]$orthFactors[j] = targetFunsList[[i]]$orthFactors[j] - orthFunDensVecNorm[,j]%*%targetFunEvals[,i]
    targetFunEvals[,i] = targetFunEvals[,i] + targetFunsList[[i]]$orthFactors[j]*orthFunEvals[,j]
  }
  for(j2 in seq_len(i-1)){
    targetFunsList[[i]]$orthFactors[seq_len(j2+orthLength)] = targetFunsList[[i]]$orthFactors[seq_len(j2+orthLength)] - c(getInnerScoreRatio(targetFunEval = targetFunEvals[,i], orthFunDensVecNorm = targetFunEvals[,i-j2]*densVec/c(densVec %*% targetFunEvals[,i-j2]^2),...))
  targetFunEvals[,i] = targetFunEvals[,i] + targetFunsList[[i]]$orthFactors[orthLength+j2]*targetFunEvals[,i-j2]
  }
}
  if(normalize) {
    targetFunEvals = rowMultiply(targetFunEvals, 1/sqrt(densVec %*% targetFunEvals^2))
  }
  return(list(targetFunsList = targetFunsList, targetFunEvals = targetFunEvals))
}
