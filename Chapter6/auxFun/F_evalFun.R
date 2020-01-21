#' An evaluation function
evalFun = function(test, seqMat, B, verifDat = NULL, sigLevel = 0.1,
                   samr = TRUE, zValues = TRUE, numSignif = TRUE,...){
testRes = testFun(seqMat = seqMat, B = B, test = test, samr =samr,
                  zValues =zValues, ...)
  if(!is.null(verifDat)){
testResVerif = testFun(seqMat = verifDat, B = B, test = test, zValues = zValues, ...)
verifList = makeList(testResVerif, samr =samr)
evalList = lapply(names(verifList), function(name){
  sapply(makeList(testRes, samr = samr), evaluatePerformance, method = "none", sigLevel = sigLevel, idDA = verifList[[name]]<sigLevel)
})
names(evalList) = names(verifList)
  } else {
  evalList = sapply(makeList(testRes, samr = samr), evaluatePerformance,
                    method = "none", sigLevel = sigLevel,
                    idDA = grepl("-TP", names(testRes$pvals)),
                    numSignif = numSignif)
  }
  return(evalList)#, zvalObs = testRes$BMA$zValObs, zValsPerm = testRes$BMA$zValsMat))
}
#' A testing function
testFun = function(seqMat, B, test, samr, zValues,...){
  BMAw = testDAA(Y = seqMat$seqMat, x = seqMat$x, FC = seqMat$FC, test = test,
                 B = B, weightStrat = "LHw", zValues = zValues, ...)
  cdfVal = BMAw$cdfValObs
  pVal = 2*sapply(cdfVal, function(x){min(x, 1-x)})
  names(pVal) = names(BMAw$fdr)
  zVal = qnorm(cdfVal)
  locfdrRes = getLocFdr(zVal)
  locfdrResAsym = getLocFdr(zVal, nulltype = 3)
  theorObj = do.call(what = getTheorFdr, args = BMAw)
  oracleObj = getOracleFdr(BMAw, nullID = !grepl("-TP", names(pVal)))
  #SAM
  sam = if(samr) samAux(seqMat = seqMat, B = B, test = test) else NULL
  list(BMAw = BMAw, theorObj = theorObj, pvals = pVal,
       locfdrRes = locfdrRes, locfdrResAsym = locfdrResAsym, sam = sam,
       oracleObj = oracleObj)
}
#' A function to create a list
makeList = function(testRes, samr, norm =FALSE, bma = FALSE, oracle = TRUE){
  tmpList = with(testRes, list(
    BH = p.adjust(pvals, method ="BH"),
    theorFdr = theorObj$Fdr,
    theorfdr = theorObj$fdr,
    empFdr = locfdrRes$Fdr,
    empfdr = locfdrRes$fdr,
    BMAwFdr = BMAw$Fdr,
    BMAwfdr = BMAw$fdr
  ))
  if(!is.null(testRes$locfdrResAsym)){
    tmpList$empAsymFdr = testRes$locfdrResAsym$Fdr
    tmpList$empAsymfdr = testRes$locfdrResAsym$fdr
  }

  if(bma){
    tmpList = within(tmpList, {
      BMAFdr = BMA$Fdr
      BMAfdr = BMA$fdr})
  }
  if(norm){
    tmpList = within(tmpList, {
      BMAFdrNorm = testRes$BMANorm$Fdr
      BMAfdrNorm = testRes$BMANorm$fdr
      BMAwFdrNorm = testRes$BMAwNorm$Fdr
      BMAwfdrNorm = testRes$BMAwNorm$fdr
    })
  }
  if(oracle){
    tmpList = within(tmpList, {
      oracleFdr = testRes$oracleObj$Fdr
      oraclefdr = testRes$oracleObj$fdr
    })
  }
  if(samr) tmpList$sam = testRes$sam
  return(tmpList)
}