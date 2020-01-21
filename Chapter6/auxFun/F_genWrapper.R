#' A function to wrap the data generation
genWrapper = function(n, p, estparamList, Template, corStrat, foldChange, p0,
                      intrStrat, rep, empirical = FALSE){
  x = as.numeric(sample(estparamList[[Template]]$x, n))-1
  p = min(p,length(estparamList[[Template]]$coef))#Adapt sample sizes to availability
  idTax = order(estparamList[[Template]]$coef, decreasing = TRUE)[seq_len(p)] #Use most abundant taxa
  idSam = sample(seq_along(estparamList[[Template]]$Libs), n)
  NBcoef = estparamList[[Template]]$coef[idTax]
  NBcoef = NBcoef - log(sum(exp(NBcoef)))
  NBtheta = estparamList[[Template]]$theta[idTax]
  Libs = estparamList[[Template]]$Libs[idSam]
  FC = if(grepl("FC", corStrat)) estparamList[[Template]]$FC[idSam] else 1
  Sigma = if(grepl("Cor", corStrat)) estparamList[[Template]]$Sigma else NULL
  idDA = sample(names(NBcoef), round((1-p0)*p)) #Sample IDs of differential abundance taxa
  intFC = introduceFC(NBcoef = NBcoef, FC = FC, idDA = idDA, foldChange = foldChange, intrStrat = intrStrat, xLog = x==1) #Introduce DA
  FC = intFC$FC
  coefMat = matrix(NBcoef, n, p, byrow = TRUE)
  coefMat[x==1,] = matrix(intFC$NBcoef, sum(x==1), p, byrow = TRUE)
  seqMat = genNB(p = p, mu = exp(coefMat)*Libs, NBtheta = NBtheta, x = NULL,
                 libSizes = NULL, beta = NULL, Sigma = Sigma, idDA = NULL,
                 n = n, taxNames = names(intFC$NBcoef), samNames = x,
                 empirical = empirical)
  list(seqMat = seqMat, FC = FC, x = x)
}
#' Plasmodes
genWrapperSimSeq = function(physeq, samSize, groupName, concName, lfdr, p0, p, Cor){
  if(!taxa_are_rows(physeq)){
    physeq=t(physeq)
  }
  p = min(p, ntaxa(physeq))
  counts = physeq@otu_table@.Data
  taxID = order(rowSums(counts), decreasing = TRUE)[seq_len(p)]
  counts = counts[taxID,]
  treatment = get_variable(physeq, groupName)
  idNA = is.na(treatment)
  treatment = treatment[!idNA]
  counts = counts[, !idNA]
  plasmode = SimData(counts, treatment, sort.method = "unpaired", k.ind = samSize, n.diff = round(p*(1-p0)), weights = lfdr[taxID],
                     n.genes = p, norm.factors = colSums(counts), switch.trt = table(treatment)[1]<table(treatment)[2])
  rownames(plasmode$counts)[plasmode$DE.ind] = paste0(rownames(plasmode$counts)[plasmode$DE.ind], "-TP")
  list(seqMat = t(plasmode$counts), FC = if(Cor=="Cor") 1 else sample(get_variable(physeq, concName), samSize*2), x = plasmode$treatment)
}
#' The mock data
genWrapperMock = function(n, p, TemplatePhy, rep, FCname, Cor){
  idSam = sample(sample_names(TemplatePhy), n)
  TemplatePhy = prune_taxa(TemplatePhy, taxa = taxa_names(TemplatePhy)[order(taxa_sums(TemplatePhy))[seq_len(p)]])
  TemplatePhy =  prune_samples(TemplatePhy, samples = idSam)
  list(seqMat = otu_table(TemplatePhy)@.Data, FC = if(Cor=="Cor") 1 else get_variable(TemplatePhy, varName = FCname), x = sample(rep(c(0,1), each = n/2)))
}
#' Evalutation-verification
genWrapperEV = function(n, p, TemplatePhy, rep, FCname, groupName, Cor, largerFactor = 4){
  #Most abundant taxa
  TemplatePhy = prune_taxa(TemplatePhy, taxa = taxa_names(TemplatePhy)[order(taxa_sums(TemplatePhy))[seq_len(p)]])
  x = get_variable(TemplatePhy, varName = groupName)
  idSamLog = x==x[1]
  idSam0 = sample(sample_names(TemplatePhy)[idSamLog], n/2*largerFactor)
  idSam1 = sample(sample_names(TemplatePhy)[!idSamLog], n/2*largerFactor)
  #Sample evenly per treatment group
  idEval = c(sample(idSam0, n/2), sample(idSam1, n/2))
  idVerif = c(idSam0, idSam1)[!c(idSam0, idSam1) %in% idEval]
  TemplatePhyEval = prune_samples(TemplatePhy, samples = idEval)
  TemplatePhyVerif = prune_samples(TemplatePhy, samples = idVerif)
  list(
    eval = list(seqMat = otu_table(TemplatePhyEval)@.Data, FC = if(Cor=="Cor") 1 else get_variable(TemplatePhyEval, varName = FCname), x = as.integer(get_variable(TemplatePhyEval, varName = groupName))-1),
    verif = list(seqMat = otu_table(TemplatePhyVerif)@.Data, FC = if(Cor=="Cor") 1 else get_variable(TemplatePhyVerif, varName = FCname), x = as.integer(get_variable(TemplatePhyVerif, varName = groupName))-1)
    )
}
