#' An auxiliary function to estimate NB parameters from a data matrix or phyloseq object
estNBparams = function(Y, X, prevCutOff = 1){
  Y = Y[rowSums(Y)>0, colMeans(Y==0)<prevCutOff]
Libs =  rowSums(Y)
tmpFit = apply(Y, 2, function(y){
  nbFit = glm.nb2(y = y, reg = X, s = Libs)
  list(theta = nbFit$theta, coef = nbFit$betas)
})
coef = sapply(tmpFit, function(x){x$coef})
list(Libs = Libs, theta = sapply(tmpFit, function(x){x$theta})[!is.na(coef)],
     coef = coef[!is.na(coef)], X = X)
}
estNBparamsPhylo = function(physeq, phyVars, round = FALSE,...){
  otuTab = if(taxa_are_rows(physeq)) t(otu_table(physeq)@.Data) else otu_table(physeq)@.Data
  if(round) otuTab = round(otuTab)
  modelMat = model.matrix(data = data.frame(sample_data(physeq)), object = formula(paste("~",paste(phyVars, collapse ="+"))))
  estNBparams(Y = otuTab[rownames(modelMat),], X = modelMat,...)
}
