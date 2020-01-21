#' A function to plot the results of the simulation studies using boxplots.
#'
resPlot = function(resultsList, varList, allowedList, rows, cols, y = "value",
                   x = "Multiplicity", colour = "Multiplicity",
  name = if("Criterion" %in% varList) allowedList[varList=="Criterion"] else "",
  showMeans = TRUE, sigLevel = 0.1,
  yIntercept = ifelse("FDP" %in% allowedList, 0.1,0.5), showFdr = FALSE, palette  = "Paired",
  legendName = "Multiplicity correction",...){
  #Filter the data
  data = subsetDf(resultsList, varList = varList, allowedList = allowedList)
  data = droplevels(data[data$Criterion != "pos",])
  Plot = ggplot(data = data, mapping = aes_string(x = x, y = y, colour = x)) +
    geom_boxplot() +
    facet_grid(reformulate(cols, rows)) +
    scale_y_continuous(name = name) +
    theme_bw() +
    theme(axis.text.x = element_blank(), axis.ticks = element_blank()) +
    scale_colour_brewer(palette = palette, name = legendName)
    if(rows=="Criterion") {
    Plot = Plot + geom_hline(show.legend = FALSE, data = data.frame(h = c(0.5, sigLevel), Criterion = c("Sensitivity","FDP"  )), mapping = aes(yintercept = h), col ="darkgreen", linetype ="dashed", size = 0.5)
    } else {
    Plot = Plot + geom_hline(yintercept = yIntercept, col ="darkgreen", linetype ="dashed", size = 0.5)
    }
  Form = formula(paste(y, "~", paste(collapse = "+", c(x, rows, cols))))
    if(showMeans) {
      #Calculate means
      aggData = aggregate(data = data, Form, FUN = mean)
      Plot = Plot + geom_point(data = aggData, shape = 23, col = "black")
    }
  if(showFdr){
    allowedList = lapply(allowedList, function(x){
      if(x=="FDP") return(list("FDP", "pos")) else return(x)})
    dataW = subsetDf(resultsList, varList = varList,
                     allowedList = allowedList)
    w = getFdrWeights(dataW, c(x, rows, cols))
    Plot = Plot + geom_point(data = w, shape = 23, col = "blue") + xlab("")
  }
  Plot
}