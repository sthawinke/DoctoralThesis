#' A function to add the empirical false discovery rate to lists lacking it
addLocFdr = function(evalRes){
  locFdr = sapply(getLocFdr(evalRes$zvalObs), evaluatePerformance, method = "none", sigLevel = 0.1, idDA = grepl("-TP", names(evalRes$zvalObs)))
  list(evalList = cbind(evalRes$evalList, locFdr), zvalObs = evalRes$zvalObs)
}