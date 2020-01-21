#' Plot pseudoF statistics
plotPseudoF = function(samClusList, yVar = "pseudoF", xVar = "Method",
                       cols = "Template", size = 0.5, meanSize = 1.2,
                       rows = NULL, scales = "fixed",
                       colNames = c(yVar, xVar, "rep", cols, rows)){
    samAmolt = melt(samClusList)
    names(samAmolt) = colNames
    samAmolt$Method = factor(samAmolt$Method, labels = methodsLabels,
                             levels = methodsLevels, ordered = TRUE)
    Plot = ggplot(samAmolt, aes_string(y = yVar, x = xVar, col = xVar)) +
        geom_boxplot(size = size) +
        theme_bw() +
        theme(axis.text.x = element_blank(), axis.ticks = element_blank()) +
        facet_grid(reformulate(cols, rows), scales = scales) +
        ylab("Pseudo F-statistic")
    #Means
    Form = formula(paste(yVar, "~", paste(collapse = "+", c(xVar, cols, rows))))
    aggData = aggregate(data = samAmolt, Form, FUN = mean)
    Plot = Plot + geom_point(data = aggData, shape = 23, col = "black",
                             size = meanSize)
    #Colour
    Plot = Plot + scale_colour_manual(values = methodsColors[methodsLabels %in%
                                                levels(droplevels(samAmolt$Method))])
    Plot
}