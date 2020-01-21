#' A function to generate correlated NB data, given a covariance matrix
#' @param n number of observations
#' @param mu means of NB distribution
#' @param Sigma a positive definite covariance or correlation matrix
#' @param ks overdispersion parameters (size)
Rmvnegbin = function (n, mu, Sigma, ks, empirical, ...)
{
  require(MASS)
  Cor <- cov2cor(Sigma)
  if (missing(mu))
    stop("mu is required")
  if (dim(mu)[2] != dim(Sigma)[2])
    stop("Sigma and mu dimensions don't match")
  if (missing(ks)) {
    SDs <- sqrt(diag(Sigma))
    ks <- unlist(lapply(1:length(SDs), function(i) .negbin_getK(mu[i],
                                                                SDs[i])))
  }
  d <- dim(mu)[2]
  normd <- MASS::mvrnorm(n, mu = rep(0, d), Sigma = Cor, empirical = empirical) #The normal-to-anything framework
  unif <- pnorm(normd)
  data <- t(qnbinom(t(unif), mu = t(mu), size = ks, ...))
  data <- .fixInf(data)
  return(data)
}

##An auxiliary function
.fixInf <- function(data) {
  # hacky way of replacing infinite values with the col max + 1
  if (any(is.infinite(data))) {
    data <-  apply(data, 2, function(x) {
      if (any(is.infinite(x))) {
        x[ind<-which(is.infinite(x))] <- NA
        x[ind] <- max(x, na.rm=TRUE)+1
      }
      x
    })
  }
  data
}
