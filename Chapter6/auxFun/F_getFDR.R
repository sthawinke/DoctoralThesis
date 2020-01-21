#' Functions to extract the fdr and drop unused levels
getFDR = function(df, vars){
  formula = formula(paste("value ~", paste(vars, collapse = "+")))
  aggregate(data = df, formula, FUN = mean)
}
filterLevelsPlot = function(df, levels, labels){
  df = df[df$Multiplicity %in% levels,]
  df$Multiplicity = factor(df$Multiplicity, levels = levels,
                           labels = labels, ordered = TRUE)
  df
}