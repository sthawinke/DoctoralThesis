#' A subsetting function for plotting
subsetDf = function(df, varList = NULL, allowedList = NULL){
  if(is.null(varList)) df
  if(length(varList)!=length(allowedList)) stop("Number of variables must match with number of levels!\n")
  for (i in seq_along(varList)){
    df = df[df[[varList[[i]]]] %in% allowedList[[i]],]
  }
  return(df)
}