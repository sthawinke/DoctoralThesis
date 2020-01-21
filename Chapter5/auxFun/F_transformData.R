#' Transform the data as required
transformData  =function(datList, what){
    mapply(datList, what, SIMPLIFY = FALSE,
           FUN = function(dat, wat){
        if(wat == "pca"){
            scale(dat, center = TRUE, scale = FALSE)
        } else if(wat == "ca"){
            rs = rowSums(dat); cs =  colSums(dat)
            E = outer(rs,cs)/sum(dat)
            (dat-E)/sqrt(E)
        } else if(wat == "clr"){
            clrTransformMat(dat)
        } else {stop("Option unknown!\n")}
    })
}