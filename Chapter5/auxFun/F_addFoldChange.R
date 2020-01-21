#' Introduce fold change in relative abundances composition
addFoldChange = function(rhos, fc, H1frac = TPR, compensate=FALSE)
{
    if(fc==1) {return(rhos)}
    nTaxa = length(rhos)
    if(compensate){
        nOTUsUp = round(nTaxa*H1frac*(1/(fc+1))) #Upregulated taxa
        nOTUsDown = round(nTaxa*H1frac-nOTUsUp) #Downregulated taxa
        someNegatives=TRUE; rhosOrig  = rhos
        #Make sure all abundances are positive
        while(someNegatives){
        rhos =rhosOrig
        OTUids = sample(names(rhos), nOTUsUp + nOTUsDown, replace=FALSE)
        OTUidUps = OTUids[1:nOTUsUp]
        OTUidDowns = OTUids[(nOTUsUp+1):(nOTUsDown+nOTUsUp)]
        rhos[OTUidUps] = rhos[OTUidUps]*fc # Add fold change up
        rhos[OTUidDowns] = rhos[OTUidDowns]*(1-sum(rhos[OTUidUps])-sum(rhos[!(names(rhos) %in% OTUids)]))/sum(rhos[OTUidDowns]) #And compensate the downs. This way the average FC is 5 in both directions and the TN taxa are really left untouched
        indTPup <- names(rhos) %in% OTUidUps
        newTaxaNamesUp <- paste0(names(rhos)[indTPup], "-TPup")
        indTPdown <- names(rhos) %in% OTUidDowns
        newTaxaNamesDown <- paste0(names(rhos)[indTPdown], "-TPdown")
        names(rhos)[indTPup] <- newTaxaNamesUp
        names(rhos)[indTPdown] <- newTaxaNamesDown
        someNegatives = !all(rhos > 0)
        }
    } else {
        nOTUs = round(nTaxa*H1frac) #DA taxa
        OTUids = sample(names(rhos), nOTUs, replace=FALSE)
        rhos[OTUids] = rhos[OTUids]*fc # Add fold change up
        indTP <- names(rhos) %in% OTUids
        newTaxaNames <- paste0(names(rhos)[indTP], "-TPup")
        names(rhos)[indTP] <- newTaxaNames
    }
    rhos/sum(rhos) #Renormalize.
}