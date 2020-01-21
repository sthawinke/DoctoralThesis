#' A function to extract deviance residuals
getDeviances = function(params, Y){
  Y = Y[, names(params$theta)]
  mu = params$Libs * exp(params$X %*% params$coef)
  thetaMat = matrix(params$theta, nrow(Y), ncol(Y), byrow =TRUE)
  tmpMat = sqrt(2*(Y*log(Y/mu)-(Y+thetaMat)*log((1+Y/thetaMat)/(1+mu/thetaMat))))*sign(Y-mu)
  tmpMat[Y==0] = -sqrt((2*thetaMat*log(1+mu/thetaMat))[Y==0]) #zero observations are always smaller than the mean
  tmpMat
}