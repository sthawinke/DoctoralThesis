#' Functions related to construction of plasmodes with SimSeq
getLFDRs = function(lfdrList, test, template, methodsToUse = c("theorfdr", "empfdr","BMAfdr","BMAwfdr")){
  #It somehow feels like we need geometric means here, since we deal with probabilities
  TDR = 1-apply(Reduce(f = cbind, lfdrList[[test]][[template]][methodsToUse]), 1, function(x){
    x=x[!is.na(x)]
    exp(mean(log(x)))
  })
  TDR
}