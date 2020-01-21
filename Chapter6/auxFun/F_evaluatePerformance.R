#' A function to evaluate performance
evaluatePerformance = function(pVal, idDA = character(0), sigLevel = 0.05,
                               method ="BH", numSignif = FALSE){
pAdj = p.adjust(pVal,method = method)
pAdj[is.na(pAdj)] = 1
if(is.null(names(pVal))) {names(pVal)= seq_along(pVal)}
pos = names(pVal)[pAdj < sigLevel]
neg = names(pVal)[!names(pVal) %in% pos]
if(is.logical(idDA)){idDA = names(pVal)[idDA]}
truePos = sum(pos %in% idDA)
falsePos = length(pos) - truePos
trueNeg = sum(!neg %in% idDA)
falseNeg = length(neg) - trueNeg
Sens = truePos/(truePos+falseNeg)
# Spec = trueNeg/(trueNeg+falsePos)
FDP = if(length(pos)) falsePos/(falsePos + truePos) else 0
if(numSignif){
  return(c(Sens = Sens, FDP = FDP, pos = length(pos)))
  } else {
return(c(Sens = Sens, FDP = FDP))
  }
}