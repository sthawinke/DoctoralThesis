#' An auxiliary function to estimate NB parameters from a data matrix or phyloseq object
estparams = function(Y, x, FC, prevCutOff = 1, minLibSize = 1e2L){
  samID = rowSums(Y)>minLibSize
Y = Y[samID, colMeans(Y==0)<prevCutOff]
Libs = rowSums(Y)
tmpFit = apply(Y, 2, function(y){
  nbFit = glm.nb2(y = y, reg = matrix(1, nrow(Y),1), s = Libs)
  list(theta = nbFit$theta, coef = nbFit$betas)
})
coef = sapply(tmpFit, function(x){x$coef})
notNA = !is.na(coef)
list(Libs = Libs,
     theta = sapply(tmpFit, function(x){x$theta})[notNA],
     coef = coef[notNA],
     x = x[samID],
    FC = FC[samID]
     )
}
estparamsPhylo = function(physeq, groupName, FCname, ...){
  otuTab = if(taxa_are_rows(physeq)) t(otu_table(physeq)@.Data) else otu_table(physeq)@.Data
  estparams(Y = otuTab, FC = get_variable(physeq, FCname), x = factor(get_variable(physeq, groupName)),...)
}
