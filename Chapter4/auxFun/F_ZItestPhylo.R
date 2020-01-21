#' A function implementing the lkelihood ratio test to compare the fit of the NB and ZINB distributions
ZItestPhylo = function(physeq, groupVar){
  physeq = prune_taxa(physeq, taxa = taxa_sums(physeq>0))
  otuTab = if(taxa_are_rows(physeq)) t(otu_table(physeq)@.Data) else otu_table(physeq)@.Data
  modelMat = model.matrix(data = data.frame(sample_data(physeq)), object = formula(paste("~", paste0(groupVar, collapse = "+"))))
  apply(otuTab[rownames(modelMat),], 2, LRTest, s = rowSums(otuTab[rownames(modelMat),]), x = modelMat)
}