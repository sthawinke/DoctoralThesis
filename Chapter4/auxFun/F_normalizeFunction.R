#' A function to normalize functions
#'
#' @param Fun A function
#' @param ... the arguments of the functions
#' @param truncQuant the quantile up to which to apprixmate the norm
#'
#' @return a function, normalized
normalizeFunction = function(Fun,densVec, evalTargetFun, ...){
  Fun$norm  = evalTargetFun^2 %*% densVec
  Fun
}
