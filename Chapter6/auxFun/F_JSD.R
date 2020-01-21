#' A function to calculate the Jensen-Shannon divergence between two probability distributions
JSD = function(p, q){ #The Jensen-Shannon divergence
m = 0.5 * (p + q)
sum(p * log(p / m)) + sum( q * log(q / m))
}
KL = function(p, q, weights = 1){ #The Kullback-Leibler divergence from p from q
  # p = p/sum(p); q = q/sum(q)
  sum(log(q / p) * weights)
}
KLemp = function(qEmpFun, pEmpDiff, sortedZvals, epsilon, weights = 1){
sum(
  log(pEmpDiff/(qEmpFun(sortedZvals)-qEmpFun(sortedZvals-epsilon)))
  * weights)
}