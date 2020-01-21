#' Estimate the parameters of a gaussian
estNormalParams = function(mat){
    apply(mat, 2, function(x){
        fit = lm(x~1)
        c("coef" = unname(fit$coef), "sd" = sigma(fit))
    })
}