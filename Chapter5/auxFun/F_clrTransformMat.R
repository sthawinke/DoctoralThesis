#The complete clr-transform of a matrix
clrTransformMat = function(mat, method = "CZM", output = "p-counts"){
    require(zCompositions)
    Xzero <- as.matrix(cmultRepl(X = mat, method = method,
                             output = output, suppress.print = TRUE)) #Correct zeroes in a Bayesian way
    # convert to proportions
    Xnorm = Xzero/rowSums(Xzero)
    # make our compositional dataset
    Xcrl <- t(apply(Xnorm, 1, function(x){log(x) - mean(log(x))}))
    #Column center
    Xsca = scale(Xcrl, center = TRUE, scale = FALSE)
    Xsca
}