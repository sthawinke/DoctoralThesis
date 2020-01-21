#' A function to estimate different thetas per treatment group. Below a wrapper for a matrix
estDiffThetas = function(y, xFac, reg, s, maxit = 200L, convTol = 1e-10){
  #If any subgroup all zero, model cannot be fitted
  if(any(tapply(y, xFac, sum)==0)) { return(list(betas = rep(NA, ncol(reg)), thetas = rep(NA, ncol(reg))))}
  iter = 1
  convergence = FALSE
  fitPrim = try(glm.nb2(y = y,reg = reg, s=s ))
  if(class(fitPrim)[[1]]=="try-error") { fitPrim = (list(betas = rep(0.001, ncol(reg)), thetas = rep(10, ncol(reg))))}
  thetas = rep(fitPrim$theta, length(unique(xFac)))
  betas = fitPrim$betas
  while(!convergence & (iter<= maxit)){
  mu = exp(reg %*% betas)*s
  betasOld = betas; thetasOld = thetas
  thetas = tapply(seq_along(y), xFac, function(i){theta.ml2(theta = 1e-4,y = y[i], mu = mu[i])})
  betas = try(nleqslv(x = betas, fn = ScoreNB, jac = JacNB, reg = reg , y = y, thetas = thetas[xFac], s = s)$x)
  if(class(betas)=="try-error") {betas = try(nleqslv(x = betas, fn = ScoreNB, jac = NULL, reg = reg , y = y, thetas = thetas[xFac], s = s)$x)}
  if(class(betas)=="try-error") {return(list(betas = rep(NA, ncol(reg)), thetas = rep(NA, length(unique(x)))))}
  iter = iter + 1
  convergence = (all(abs(1-betas/betasOld) < convTol)) & all(abs(thetas-thetasOld) < convTol)
  }
  if(!convergence){warning("No convergence achieved!\n")}
  return(list(betas = betas, thetas = thetas))
}
#Must supply x a factor
estDiffThetasMat = function(Y, x, prevCutOff = 0.95,...){
  if(!is.factor(x)){stop("Must supply a factor variable for separate dispersion estimation")}
  Y = Y[rowSums(Y)>0, colMeans(Y==0)<prevCutOff]
  # sapply(seq_len(ncol(Y)), function(i){cat(i, "\n");estDiffThetas(Y[,i], xFac = x, reg = model.matrix( ~x, data.frame(x=x)), s = rowSums(Y),...)})
apply(Y, 2, estDiffThetas, xFac = x, reg = model.matrix( ~x, data.frame(x=x)), s = rowSums(Y),...)
}
estDiffThetasPhylo = function(physeq, groupVar, convFactor = FALSE){
  if(length(groupVar)>1){stop("Supply only one grouping variable!\n")}
  otuTab = if(taxa_are_rows(physeq)) t(otu_table(physeq)@.Data) else otu_table(physeq)@.Data
  variab = if(convFactor) factor(get_variable(physeq, groupVar)) else get_variable(physeq, groupVar)
  otuTab = round(otuTab[!is.na(variab),])
  variab = variab[!is.na(variab)]
  estDiffThetasMat(otuTab, variab)
}