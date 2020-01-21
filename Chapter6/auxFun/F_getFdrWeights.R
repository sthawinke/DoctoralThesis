#' A function to calculate the Fdr
#' @param df A dataframe with at least variables "value" and "Criterion"
#' @param vars Conditioning variables

getFdrWeights = function(df, vars){
  formula = formula(paste("value ~", paste(vars, collapse = "+")))
  weights = aggregate(data = df[as.character(df$Criterion) == "pos", ],
                      formula, FUN = function(x){x/sum(x)})
  weights$Criterion = "FDP"
  dfFDP = df[as.character(df$Criterion) == "FDP", ]
  #If no discoveries, no weight (NA)
  Val = as.numeric(as.character(
    sapply(seq_len(nrow(weights)), function(i){
    subDf = dfFDP[apply(MARGIN = 1, FUN = all,
                     sapply(vars, function(Var){
                       as.character(dfFDP[[Var]]) ==
                         as.character(weights[[Var]])[i]})),]
    sum(subDf$value * weights$value[i,])
  })))
  weights$value = Val

  return(weights)
}