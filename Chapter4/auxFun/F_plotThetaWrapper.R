# A function to plot the different estimates
plotEstimates = function(thetasOne, thetasMat, cex = 0.8,log = "xy",...){
  int = intersect(names(thetasOne), colnames(thetasMat))
  thetasMat[thetasMat==0] = NA
  plot(thetasOne[int], thetasMat[1,int], xlab = "Single dispersion estimates", ylab = "Group-wise dispersion estimates", ylim = range(thetasMat,na.rm = TRUE),cex = cex,log = log,...)
  abline(0,1)
  legend("bottomright", col = seq_len(nrow(thetasMat)), legend = rownames(thetasMat), pch = 1)
  for (i in 2:nrow(thetasMat)) {points(thetasOne[int], thetasMat[i,int], col = i, cex = cex)}
}
# A function to extract the estimates
getThetaMat = function(Params, physeq, Var){
  x= get_variable(physeq, Var)
  tmp = sapply(unique(x), function(y){
    sapply(Params, function(z){
      z$thetas[y]
    })
  })
  colnames(tmp) = unique(x)
  rownames(tmp) = names(Params)
  tmp
}
# A wrapper for both
plotThetaWrapper = function(paramList, paramListAll, ...){
  plotEstimates(paramList$theta, paramListAll$theta, ...)
}