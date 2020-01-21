# A function to get the log-likelihood for the complex  and single dispersion case
getLogLikNB = function(y, thetas, betas, reg, s, id = rownames(reg)){
  mu = c(exp(reg %*% betas)*s)
  sum(dnbinom(y[id], mu = mu, size = thetas, log = TRUE))
}
getLogLiksParams = function(params, paramListAll, physeq, var){
  int = intersect(names(params$theta), colnames(paramListAll$theta))
  otuTab = round(if(taxa_are_rows(physeq)) t(otu_table(physeq)@.Data) else otu_table(physeq)@.Data)
  sapply(int, function(i){
    single = getLogLikNB(otuTab[,i], thetas = params$theta[i], betas = params$coef[,i], reg = params$X, s = params$Libs)
    all = getLogLikNB(otuTab[,i], thetas = paramListAll$theta[get_variable(physeq, var)[sample_names(physeq) %in% rownames(params$X)],i], betas = paramListAll$coef[,i], reg = params$X, s = params$Libs)
    c(all = all, single = single)
  })
}
#Perform the likelihood ratio test
LRT = function(llMat, df){
  pchisq(2*(llMat[1,]-llMat[2,]), df = df, lower.tail = FALSE)
}