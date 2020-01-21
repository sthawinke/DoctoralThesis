#' A function that performs a likelihood ratio test, returning a p-value
LRTest = function (y, x, s, df=1)
{
  require(pscl);require(VGAM)
  y = round(y)
  NBfit = glm.nb2(y = y, reg = x, s = s)
  NBll = sum(dnbinom(log = TRUE, x = y, size = NBfit$theta, mu = exp(x%*%NBfit$beta)*s))
  ZINBfit = try(zeroinfl(y~offset(log(s))+x-1|1,dist="negbin", y = FALSE, model = FALSE),silent=TRUE)
  ZINBll <- if(class(ZINBfit)=="try-error") NA else max(logLik(ZINBfit),sum(dzinegbin(log = TRUE, x = y, size = NBfit$theta, munb = exp(x%*%ZINBfit$coef$count)*s, pstr0 = 0)))
  pval <- stats::pchisq(2*(ZINBll-NBll), df = df, lower.tail = FALSE)
  return(pval)
}