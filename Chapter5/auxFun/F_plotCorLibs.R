#' Plot the correlations as boxplots
#' @param corLibs The calculated correations
#' @param Dims the dimensions resuired
#' @param yVar,xVar The x and y-variables
#' @param varNames the variable names
#'
#' @return The ggplot object
plotCorLibs = function(corLibs, Dims = 1:2, yVar = "Correlation", xVar = "Method",
                       varNames = c("Correlation", "Dimension", "Method", "Dataset", "Repeat"),
                       rows = "Dimension", cols = "Dataset",
                       methodsKeep = methodsClusters,
                       pointSize = 1){
    moltCor = melt(corLibs)
    names(moltCor) = varNames
    moltCor$Method = factor(moltCor$Method, labels = methodsLabels,
                            levels = methodsLevels, ordered = TRUE)
    moltCor = droplevels(moltCor[moltCor$Method %in% methodsClusters,])
    if(cols=="Dataset") {
        # moltCor$Dataset = factor(moltCor$Dataset,
        #                         levels = c(grep(unique(moltCor$Dataset),
        #                                         invert = TRUE, value = TRUE,
        #                                         pattern = "overall"),
        #                                    "overall"), ordered = TRUE)
        moltCor$Dataset = droplevels(factor(moltCor$Dataset,
                                 levels = corLevels, labels = corLabels, ordered = TRUE))
    }
    moltCor = moltCor[moltCor$Dimension %in% paste0("Dim", Dims),]
    Plot = ggplot(data = moltCor, aes_string(y = yVar, x = xVar, colour = xVar)) +
        geom_boxplot() + theme_bw() +
        facet_grid(reformulate(rows, cols)) +
        geom_hline(yintercept = 0, linetype = "dashed", colour = "grey") +
        theme(axis.text.x = element_blank(), axis.ticks = element_blank())
    #Means
    Form = formula(paste(yVar, "~", paste(collapse = "+", c(xVar, rows, cols))))
    aggData = aggregate(data = moltCor, Form, FUN = mean)
    Plot = Plot + geom_point(data = aggData, shape = 23, col = "black", size = pointSize)
    #Colour
    Plot = Plot + scale_colour_manual(values = methodsColors[methodsLabels %in%
                                                                 levels(droplevels(moltCor$Method))])
    Plot

}