#' The moment generating function for the NB, returns the first four raw moments
mgfNB = function(mu, theta){
  c(mu,
    mu^2*(1+1/theta)+mu,
    (mu^2*theta^2 + 3*mu^2*theta + 3*mu*theta^2 + 2*mu^2 + 3*mu*theta +
       theta^2)*mu/theta^2,
    (mu^3*theta^3 + 6*mu^3*theta^2 + 6*mu^2*theta^3 + 11*mu^3*theta +
    18*mu^2*theta^2 + 7*mu*theta^3 + 6*mu^3 + 12*mu^2*theta + 7*mu*theta^2 +
       theta^3)*mu/theta^3)
}
getFifthMom = function(mu, od){
  (mu^4*od^4 + 10*mu^4*od^3 + 10*mu^3*od^4 + 35*mu^4*od^2 + 60*mu^3*od^3 + 25*mu^2*od^4 + 50*mu^4*od + 110*mu^3*od^2 + 75*mu^2*od^3 + 15*mu*od^4 + 24*mu^4 + 60*mu^3*od + 50*mu^2*od^2 + 15*mu*od^3 + od^4)*mu/od^4}