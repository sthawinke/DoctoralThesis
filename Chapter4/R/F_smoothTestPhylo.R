smoothTestPhylo = function(physeq, groupVar, single = TRUE, ...){
otuTab = if(taxa_are_rows(physeq)) t(otu_table(physeq)@.Data) else otu_table(physeq)@.Data
modelMat = model.matrix(data = data.frame(sample_data(physeq)), object = formula(paste("~", paste0(groupVar, collapse = "+"))))
smoothTestMat(otuTab[rownames(modelMat),], modelMat, single = single, xFac = if(single) NULL else factor(get_variable(physeq, groupVar)),...)
}
