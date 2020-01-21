#' Extract pseudoF's
getPseudoF = function(x, y){
    # if(!identical(y[[1]]$variable, y[[2]]$variable)){
    #     stop("Variables do not match!")
    # }
    tmp = lapply(x$rows, function(z, group){
        if(is.list(z)){
            NULL #c(NaN, sapply(z, pseudoF, clusters = group))
        } else {
            pseudoF(z, clusters = group)#, NaN, NaN)
        }
        #names(tmp) = c("overall", names(y))
    }, group = y[[1]]$variable)
    tmp[!sapply(tmp, is.null)]
}