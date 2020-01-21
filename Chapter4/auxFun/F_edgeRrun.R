#' EdgeR
edgeRrun = function(otuTab, desmat, ranX){
dge <- DGEList(counts = t(otuTab), group = ranX)
dge$samples$norm.factors <- rowSums(otuTab)
dgeW <- try(estimateGLMRobustDisp(y = dge, design = desmat), silent = TRUE)
if(inherits(dgeW, "try-error")){
  pValsedgeR = rep(NA, ncol(otuTab))
} else {
  glmFit <- glmQLFit(y = dgeW, dispersion = dgeW$tagwise.dispersion,
                     robust = TRUE, design = desmat)
  glmRes <- glmQLFTest(glmFit, coef = "ranX1")
  pValsedgeR <- glmRes$table$PValue
}
names(pValsedgeR) = colnames(otuTab)
return(pValsedgeR)
}