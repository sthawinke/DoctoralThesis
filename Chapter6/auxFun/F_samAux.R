samAux = function(seqMat, B, test){
  samID = rowSums(seqMat$seqMat)>0
  seqMat$seqMat = seqMat$seqMat[samID, colSums(seqMat$seqMat)>0]
  #Construct q-value table to fit into the pipeline
  out = rep(1, ncol(seqMat$seqMat))
  names(out) = colnames(seqMat$seqMat)
  res = try(SAM(x = t(seqMat$seqMat/rowSums(seqMat$seqMat)* if(length(seqMat$FC)==1) seqMat$FC else seqMat$FC[samID]), y = seqMat$x[samID]+1, testStatistic = if(test=="t.test") "standard" else "wilcoxon", nperms = B, fdr.output = 0.1, return.x = FALSE, resp.type = "Two class unpaired"), silent = TRUE)
  if(class(res)=="try-error") {return(out)} #Failed fit => No significance declared
  if (res$siggenes.table$ngenes.up == 0 & res$siggenes.table$ngenes.lo == 0){
    out
  } else {
sigTable = rbind(res$siggenes.table$genes.lo, res$siggenes.table$genes.up)
out[as.integer(gsub("g", "",sigTable[,"Gene ID"]))] = as.double(sigTable[,"q-value(%)"])/100
  }
  out
}