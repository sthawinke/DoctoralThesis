#' A function to introduce fold changes
introduceFC = function(NBcoef, FC, idDA, foldChange, intrStrat, xLog){
  if(foldChange==1){return(list(NBcoef = NBcoef, FC = FC))}
  nTaxa = length(idDA)
  rhos = exp(NBcoef)
  if(intrStrat=="seq"){
    nOTUsUp = round(nTaxa * (1/(foldChange + 1)))  #Upregulated taxa
    nOTUsDown = nTaxa - nOTUsUp  #Downregulated taxa
    OTUidUps = idDA[1:nOTUsUp]
    OTUidDowns = idDA[(nOTUsUp + 1):(nOTUsDown + nOTUsUp)]
    while(sum(rhos[OTUidUps] * (foldChange-1))> sum( rhos[OTUidDowns])){ #Avoid negative abundances
  idDA = sample(idDA)
  OTUidUps = idDA[1:nOTUsUp]
  OTUidDowns = idDA[(nOTUsUp + 1):(nOTUsDown + nOTUsUp)]
    }
    rhos[OTUidUps] = rhos[OTUidUps] * foldChange  # Add fold change up
    rhos[OTUidDowns] = rhos[OTUidDowns] * (1 - sum(rhos[OTUidUps]) - sum(rhos[!(names(rhos) %in%
                       c(OTUidUps, OTUidDowns))]))/sum(rhos[OTUidDowns]) #And compensate the downs. This way the average FC is the same in both directions and the TN taxa are really left untouched
    indTPup = names(rhos) %in% OTUidUps
    newTaxaNamesUp = paste0(names(rhos)[indTPup], "-TPup")
    indTPdown = names(rhos) %in% OTUidDowns
    newTaxaNamesDown = paste0(names(rhos)[indTPdown], "-TPdown")
    names(rhos)[indTPup] = newTaxaNamesUp
    names(rhos)[indTPdown] = newTaxaNamesDown
    NBcoef = log(rhos/sum(rhos)) #TESTED
    } else if (intrStrat=="both"){
    rhos[idDA] = rhos[idDA]*foldChange
    FC[xLog] = FC[xLog]*sum(rhos)
    names(rhos)[names(rhos) %in% idDA] = paste0(names(rhos)[names(rhos) %in% idDA], "-TP")
    NBcoef = log(rhos/sum(rhos))
    }
  return(list(NBcoef = NBcoef, FC = FC))
}