#' Plot feature correlations
plotFeatCor = function(featCor, xVar = "Method", yVar = "zValue",
                       cols = "Template", rows = "which"){
    featCor = lapply(featCor,
                     function(x){
                         lapply(x, function(y){
                            rbind(sapply(y, identity),
                                  "reference" = rnorm(length(y)))
                             })
                         })
    #Add a reference
    moltFeat = melt(featCor)
    names(moltFeat) = c(xVar, rows, yVar, "replicate", cols)
    moltFeat$Method = factor(moltFeat$Method, labels = methodsLabels,
                            levels = methodsLevels, ordered = TRUE)
    moltFeat$which = factor(moltFeat$which, levels = c("between", "all"),
                            labels = c("Between views", "All taxa"), ordered = TRUE)
    Plot = ggplot(data = moltFeat, aes_string(x = xVar, y  = yVar, col = xVar)) +
        geom_boxplot() +theme_bw() +
        facet_grid(reformulate(cols,rows)) +
        geom_hline(yintercept = 0, linetype = "dashed", colour = "grey") +
        theme(axis.text.x = element_blank(), axis.ticks = element_blank()) +
        ylab("Wilcoxon rank sum test statistic")
    #Means
    Form = formula(paste(yVar, "~", paste(collapse = "+", c(xVar, rows, cols))))
    aggData = aggregate(data = moltFeat, Form, FUN = mean)
    Plot = Plot + geom_point(data = aggData, shape = 23, col = "black", size = 1)
    #Colour
    Plot = Plot + scale_colour_manual(values = methodsColors[methodsLabels %in%
                                                                 levels(droplevels(moltFeat$Method))])
    Plot
}