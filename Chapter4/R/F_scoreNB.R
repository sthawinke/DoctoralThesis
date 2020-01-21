ScoreNB = function(betas, reg, y, thetas, s){
  mu = c(exp(reg %*% betas)*s)
  crossprod(reg,((y-mu)/(1+mu/thetas)))
}
JacNB = function(betas, reg, y, thetas, s){
  mu = c(exp(reg %*% betas)*s)
  -crossprod(reg*c((1+y/thetas)*mu/(1+mu/thetas)^2), reg)
}