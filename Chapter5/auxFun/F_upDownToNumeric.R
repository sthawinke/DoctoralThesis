#' Convert a matrix with "up" and "down" entries to 1 and -1
#' @param mat The matrix
#' @return the converted matrix
upDownToNumeric = function(mat){
    upId = mat == "up"
    outMat = matrix(-1, nrow(mat), ncol(mat), dimnames = dimnames(mat))
    outMat[upId] = 1
    outMat
}