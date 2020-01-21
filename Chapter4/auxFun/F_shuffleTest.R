# Assign samples to random groups and perform significance testing
shuffleTest = function(physeq, variable, nRep, nCores = 1, tag){
  physeq = if(taxa_are_rows(physeq)) t(physeq) else physeq
  physeq = prune_samples(x = physeq,
                         samples = !is.na(get_variable(physeq, variable)))
  physeq = prune_taxa(physeq, taxa = taxa_sums(physeq)>0)
  otuTab = as(otu_table(physeq), "matrix")
  otuTab = otuTab[, colSums(otuTab)>=1]
  logLibSizes = log(rowSums(otuTab)) #Library sizes as offsets
  realX = get_variable(physeq, variable) #The known grouping variable
  resList = mclapply(mc.cores = nCores, mc.preschedule = FALSE, seq_len(nRep), function(i){
    smallFile = paste0("shuffleRes/", tag, "/pValsList", i, ".RData")
    if(!file.exists(smallFile)){
    ranX = factor(sample(c(0,1), replace = TRUE, nsamples(physeq))) #The random variable
    #nb ML
    nbFits = lapply(seq_len(ncol(otuTab)), function(j){
      try(glm.nb(otuTab[,j] ~ realX + ranX + offset(logLibSizes)), silent = TRUE)
    })
    names(nbFits) = colnames(otuTab)
    nbFits = nbFits[!sapply(nbFits, inherits, "try-error")]
    pValsnbML = sapply(nbFits, function(fit){
      summary(fit)$coef[grep(value = TRUE, rownames(summary(fit)$coef), pattern = "ranX"), "Pr(>|z|)"]
    })
    desmat = model.matrix(~ranX + realX)
    #edgeR
    pValsedgeR = edgeRrun(otuTab, desmat, ranX)
    #DESeq2
    pValsDESeq2 = DESeq2Run(physeq, realX, ranX)
    #Friedman
    pValsFriedman = apply(otuTab/rowSums(otuTab), 2, function(y){
      prentice.test(y = y, groups =  ranX, blocks =  realX)$p.value
    })
    pValsList = cbind("nbML" = pValsnbML, "edgeR" = pValsedgeR[names(pValsnbML)],
                     "DESeq2" = pValsDESeq2[names(pValsnbML)],
                     "Prentice" =  pValsFriedman[names(pValsnbML)])
    save(pValsList, file = smallFile)
    } else load(smallFile)
    return(pValsList)
  })
  names(resList) = seq_len(nRep)
  resList
}