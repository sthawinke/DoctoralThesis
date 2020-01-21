#' Extract row scores and column loadings
#' @param DIres the compInt object
#' @param othersRes The list of other data integration objects
#'
#' @return A list with row and cols lists
extractScores = function(DIres, othersRes){
    mapply(DIres, othersRes, FUN = function(DI, others){
    rows = c(list("compInt" = if(is.null(DI))
        matrix(NA, nrow = nrow(others[["pca"]]$sampleScores),
               ncol = ncol(others[["pca"]]$sampleScores)) else DI$latentVars),
             lapply(others, function(x){x$sampleScores}))
    cols = c(list("compInt" = if(is.null(DI)) lapply(others[["pca"]]$featureLoadings, function(ii){
        matrix(nrow = nrow(ii), ncol = ncol(ii))
    }) else lapply(DI$paramEsts, function(x){t(x)[,1:2]})),
             lapply(others, function(x){x$featureLoadings}))
    # for(i in seq_along(cols)){
    #     if(is)
    #     colnames(rows[[i]]) = paste0("Dim", seq_len(ncol(rows[[i]])))
    #     for (j in seq_along(cols[[i]])){
    #         colnames(cols[[i]][[j]]) = paste0("Dim", seq_len(ncol(cols[[i]][[j]])))
    #     }
    # }
    list("rows" = rows, "cols" = cols[!sapply(cols, is.null)])
    }, SIMPLIFY = FALSE)
}