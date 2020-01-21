# SimSeq tests
testSimSeq = function(mat, variable, countMat = TRUE, ...){
    require(fdrtool)
    variable = factor(variable)
    if(countMat){
        libSizes = rowSums(mat)
    }
    pvals = apply(mat, 2, function(x){kruskal.test(x, g = variable)$p.value})
    overSum = sum(mat)
    result = apply(mat, 2, function(x){
        # avRank = length(x)/2 + 0.5
        # rankX = rank(x)
        relAbs = if(countMat) sum(x)/overSum else mean(x)
        sapply(levels(variable), function(y){
            tmp = if(countMat) sum(x[variable==y])/sum(libSizes[variable==y]) else
                mean(x[variable==y])
            ifelse(tmp > relAbs, "up", "down")
        })
        }) #Was taxon up or down regulated?
    lfdr = if(length(pvals) >= 200)
        fdrtool(pvals, statistic = "pvalue", plot = FALSE,
                verbose = FALSE)$lfdr else
            p.adjust(pvals, method = "BH")
    list(lfdr = lfdr, variable = variable, result = result,
         mat = mat, countMat = countMat)
}
testSimSeqList = function(datList, variables, countMats){
    mapply(datList, variables, countMats, FUN = testSimSeq, ...)
}