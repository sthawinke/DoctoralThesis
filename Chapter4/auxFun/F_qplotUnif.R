#' A function to make a qqplot versus the uniform distribution
qplotUnif = function(pp, ...){
  plot(sort(pp[!is.na(pp)]), seq(0,1, length.out = sum(!is.na(pp))), ..., ylab = "Expected quantiles", xlab ="Observed p-values")
  abline(0,1)
}