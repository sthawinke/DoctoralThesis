# DESseq2
DESeq2Run = function(physeq, realX, ranX){
  sample_data(physeq)$realX =realX; sample_data(physeq)$ranX =ranX
  dds <- phyloseq_to_deseq2(physeq, formula("~ realX + ranX"))
  sizeFactors(dds) <- sample_sums(physeq)
  fitMethods <- c("local", "parametric", "mean")
  for (fitRun in seq_along(fitMethods))
  {
    suppressWarnings(
      ddsDisp <- try(
        estimateDispersions(dds, fitType = fitMethods[fitRun], quiet = TRUE),
        silent = TRUE))
    if (!inherits(ddsDisp, "try-error"))
    {
      break
    } else {}
  }
  ## if no fit was successful, raise error
  if (inherits(ddsDisp, "try-error"))
  {
    pValsDESeq2 = rep(NA, ntaxa(physeq))
  } else {
  ddsRes <- nbinomWaldTest(ddsDisp)
  ddsRes <- results(ddsRes, name = grep(resultsNames(ddsRes),
                                        pattern = "ranX", value = TRUE))
  pValsDESeq2 = ddsRes[, "pvalue"]
  }
  names(pValsDESeq2) = taxa_names(physeq)
  pValsDESeq2
}