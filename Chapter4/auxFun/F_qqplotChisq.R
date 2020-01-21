qqplotChisq = function(x, df = 2,...){
  qqplot(x, qchisq(seq_along(x)/length(x), df = df), ylab = "Theoretical quantiles", xlab = "Observed quantiles", ...)
  abline(0,1)
}
