#' A function to extract tail area and local fdr from the empirical fit
getLocFdr = function(zVals, ...){
  locfdrRes = try(locfdr(zVals, plot = 0, ...), silent = TRUE)
  if(class(locfdrRes)=="try-error") {return(list(fdr = rep(1, length(zVals)),
                                                 Fdr = rep(1, length(zVals)),
                                                 failed = TRUE))}
  #Failed fit => No discoveries
  idBin = sapply(zVals, function(x){which.min(abs(x-locfdrRes$mat[, "x"]))})
  locFdr = apply(locfdrRes$mat[, c("Fdrright","Fdrleft")], 1,min)[idBin]
  names(locfdrRes$fdr) = names(locFdr) = names(zVals)
  return(list(fdr = locfdrRes$fdr, Fdr = locFdr, failed = FALSE))
}