#' generate simseq data
#'
#' @param testres The testing results
#' @param NFeat Number of features
#' @param fracFeats Fraction of DA features
#' @param samIDee sample names from previous dataset
genSimSeq = function(testres, NFeat = ncol(testres$mat), fracFeats = 0.2,
                     samIDee = NULL){
    with(testres,
         {weights = 1-lfdr
         weights[weights==0] = .Machine$double.eps
    weights = weights/sum(weights) #Normalize sampling weights
    normFactors = rowSums(mat)
    DEtaxa = sample(colnames(mat), round(NFeat*fracFeats), prob = weights)
    # The DE taxa are sampled from all taxa with the 1-lfdr as sampling weights.
    #The NDE taxa are simply the remaining ones
    whichTrt = names(sort(table(variable), decreasing = TRUE)[1])
    # The most frequent treatment group, frow which we will sample the EE taxa
    trtLeft = unique(variable)[unique(variable) != whichTrt]
    #The remaining treatment groups
    if(is.null(samIDee)) samIDee = sample(which(variable == whichTrt),
                                          nrow(mat), replace = TRUE)
    #NDA samples are taken from the largest group
    dataSimSeq = mat[samIDee,] #EE taxa
    for (i in trtLeft){
        rowInd = which(variable==i)
        dataSimSeq[rowInd, DEtaxa] = round(mat[rowInd, DEtaxa]*
                            normFactors[samIDee][rowInd]/
                                normFactors[rowInd])
        #DA taxa are taken from the other groups, adapting the sequencing depth to the most frequent one. samIDee refers to the sample indices of the new matrix
    }
    keepID = colSums(dataSimSeq) > 0
    list(data = dataSimSeq[, keepID], variable = variable,
         DEtaxa = DEtaxa[DEtaxa %in% colnames(dataSimSeq[, keepID])])}
    )
}