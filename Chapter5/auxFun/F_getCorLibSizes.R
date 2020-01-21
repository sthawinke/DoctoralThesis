#' Calculate correlation of library sizes with row scores
getCorLibSizes = function(datList, resList){
    viewNames = names(datList[[1]])
    libNames = c(viewNames, "overall")
    libSizes = lapply(datList, function(dat){
        indiLibs = lapply(dat, function(x){
            rowSums(if(is.matrix(x)) x else x$data)})
        overLib = rowSums(simplify2array(indiLibs))
        c(indiLibs, overall = list(overLib))
    })
    names(libNames) = names(libSizes) = libNames
    applyFun = function(r, l) {
        if(NROW(r)!=NROW(l)){
            a = 1
        }
        apply(r, 2, cor, l)
        }
    mapply(libSizes, resList, SIMPLIFY = FALSE,
           FUN = function(lib, res){
               lapply(libNames, function(l){
                   lapply(c(res$rows,
                            "reference" = list(matrix(rnorm(length(lib[[l]])*2),
                                ncol = 2, dimnames = list(NULL, c("Dim1", "Dim2"))))),
                         function(r){
                             if(is.list(r)) names(r) = viewNames
                                tmp = if(is.list(r)) {
                                    if(l == "overall") {
                                        rowMeans(sapply(r, applyFun, l = lib[[l]]))
                                        } else {
                                        applyFun(r[[l]], l = lib[[l]])
                                        }
                                    } else
                                    applyFun(r, l = lib[[l]])
                                as.list(tmp)
                                }
                   )
               })
           })
}