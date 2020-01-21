# A custom NB fitter
#' A function to estimate different thetas per treatment group. Below a wrapper for a matrix
glm.nb2 = function(y, reg, s, maxit = 200L, convTol = 1e-4, limit = 20L, betas = c(log(mean(y/s)),rep(1e-10, ncol(reg)-1)), theta = if(single) 0.1 else rep(0.1, length(y)), single = TRUE, xFac = NULL, se = TRUE, ...){
require(nleqslv)
  iter = 1L
  convergence = FALSE
  foo = try(silent = TRUE, while(!convergence & (iter<= maxit)){
    betasOld = betas; thetasOld = theta
    betas = nleqslv(x = betas, fn = ScoreNB, jac = JacNB, reg = reg , y = y, thetas = theta, s = s, ...)$x
    mu = exp(reg %*% betas)*s
    if(iter==1) {theta =  if(single) max(1e-7,theta.mm2(y, mu)) else tapply(seq_along(y), xFac, function(i){theta.mm2(y = y[i], mu = mu[i])})} #Start with MoM estimate
    theta = if(single) theta.ml2(theta = theta, y = y, mu = mu) else tapply(seq_along(y), xFac, function(i){theta.ml2(theta = theta[i][1],y = y[i], mu = mu[i])})[xFac]
    iter = iter + 1L
    convergence = (all(abs(betas-betasOld) < convTol)) & all(abs(theta-thetasOld) < convTol)
  })
  if(class(foo)=="try-error"){return(list(betas = rep(NA, ncol(reg)), theta =NA, vcov = matrix(NA, ncol(reg), ncol(reg))))}
  if(se) {if (single) attr(theta, "SE") <- sqrt(1/-infoOD(theta, mu, y)) else attr(theta, "SE") = tapply(seq_along(y), xFac, function(i){sqrt(1/-infoOD(theta[i], mu[i], y[i]))})}
  if(!convergence){warning("No convergence achieved!\n")}
  inf = crossprod(reg, diag(1/(1/c(mu)+1/theta))) %*% reg
  return(list(betas = betas, theta = theta, vcov = if(rcond(inf)< .Machine$double.eps) matrix(NA, ncol(reg), ncol(reg)) else solve(inf)))
}
theta.ml2 = function (y, mu, theta,...)
{
nleqslv(theta, fn = scoreOD, jac = infoOD, y = y, mu = mu, ...)$x
}
scoreOD <- function(th, mu, y) {
sum( (digamma(y+th) -digamma(th) + log(th) + 1 - log(th + mu) - (y +th)/(mu + th)))
}
infoOD <- function(th, mu, y, ...){
-sum((-trigamma(th +y) + trigamma(th) - 1/th + 2/(mu + th) - (y + th)/(mu +th)^2))
}
theta.mm2 = function(y, mu){
  length(y)/sum((y/mu - 1)^2)
}