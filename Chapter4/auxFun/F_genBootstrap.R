# A generic bootstrap function. To mitigate the problem of all zero samples for certain taxa, we repeat the sampling procedure several times (reps) and later average the p-values. If the values are still all zero, we give up on the taxon. This is a trade off between using the same parameter estimates always (as we do when boostrapping from the real data) and keeping the maximum number of taxa.
genBootStrap = function(params, reps = 5, single = TRUE, xFac = NULL, filePath = NULL, nCores = 1){
  n = length(params$Libs)
  params$coef = params$coef[,!apply(params$coef,2,anyNA)];params$theta = params$theta[!is.na(params$theta)]
  p = ncol(params$coef)
  replics = mclapply(mc.cores = nCores, mc.preschedule = FALSE, seq_len(reps), function(non){
    filPat2 = file.path(filePath, paste0("rep",non))
    if(!file.exists(paste0(filPat2, "LibsPars.RData"))){
  yBase = matrix(rnbinom(n = n*p, size = rep(params$theta, each = n), mu = exp(params$X %*% params$coef)*params$Libs), n,p, byrow = FALSE) #Crucial: only one null matrix is generated, otherwise we cannot calculate the test statistic.
  Libs = rowSums(yBase)
  names(Libs) = rownames(yBase) = names(params$Libs); colnames(yBase) = names(params$theta)
  paramEsts = apply(yBase, 2, function(y){
    nbFit = try(glm.nb2(y = y, s = Libs, reg = params$X, single = single, xFac = xFac), silent = TRUE)
    if(class(nbFit)=="try-error") {return(list(coef = rep(NA, ncol(params$X)), theta = NA))}
    return(list(coef = nbFit$betas, theta = nbFit$theta))
  })
  save(Libs, paramEsts, yBase, file = paste0(filPat2, "LibsPars.RData"))
  return(list(yBase = yBase, paramEsts = paramEsts, Libs = Libs))
    }
  })
  return(list(replics = replics, X = params$X))}
