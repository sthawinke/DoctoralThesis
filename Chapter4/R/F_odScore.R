#' A score function for the estimation of the overdispersion parameter of the negative binomial
#'
#' @param od the current value of the estimate
#' @param y a vector of counts
#' @param mu a vector of means
#' @param ... arguments for the jacobian, not used
#'
#' @return a scalar, the evaluation of the score function at od
#'
#' @details
#' The odScoreVec() version is custom for fast evaluation of a sequence going from 0 to a large integer
odScoreOld = function(y, mu,od, ySeqCum,...){
  1/od^2*log(1+mu*od) -(y+1/od)*mu/(1+mu*od) + ySeqCum[y+1]
}
odScoreVec = function(y, mu, od,...){#Function custom written for the sequence
  ySeq2 = y[-c(1:2)]-1
  tmp = c(0,0,cumsum(ySeq2/(1+ySeq2/od)))
  od^2*log(1+mu/od) -(y+od)*1/(1/mu+1/od) + tmp
}
odScoreOld2 = function(y, mu, od, ...) {
  (digamma(od +y) - digamma(od) + log(od) + 1 - log(od + mu) - (y +od)/(mu + od))
}
odScore = function(y, mu, od, ySeqCum,...) {
  (ySeqCum[y+1] + log(od) + 1 - log(od + mu) - (y +od)/(mu + od))
}
