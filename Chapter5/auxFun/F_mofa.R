# MOFA
mofaSim = function(datList, dim = 2L, maxIter = 1000L, ...){
    mofObj = createMOFAobject(lapply(datList, t))
    modOpt = getDefaultModelOptions(mofObj)
    modOpt$numFactors = dim
    trainOpt = getDefaultTrainOptions()
    trainOpt$maxiter = maxIter
    mofObjPrep = prepareMOFA(mofObj, ModelOptions = modOpt, TrainOptions = trainOpt)
    #Poisson not ideal but...
    fit = runMOFA(mofObjPrep)
    featureLoadings = lapply(getExpectations(fit, "W"), function(x){
        colnames(x) = paste0("Dim", seq_len(dim));x
    })
    sampleScores = getExpectations(fit, "Z")
    colnames(sampleScores) = colnames(featureLoadings[[1]])
    list("featureLoadings" = featureLoadings[names(datList)],
         "sampleScores" = sampleScores,
         "converged" = maxIter!= length(TrainStats(fit)$activeK))
}