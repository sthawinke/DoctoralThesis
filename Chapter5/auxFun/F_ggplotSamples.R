#' ggplot of samples from other ordination
ggplotSamples = function(samScores, variable, title, samVar = "",
                         samShape = NULL, shapeName = NULL,
                         shapeValues = 21:25, pointSize = 3){
    df = data.frame(x = samScores[,1], y = samScores[,2], col = variable)
    ggplot(aes(x = x, fill = col, y = y, shape = samShape), data = df) +
        geom_point(size = pointSize) + coord_fixed() + ggtitle(title) + theme_bw() +Theme +
        xlab("Dim1") + ylab("Dim2")  +
        (if(is.factor(variable)) scale_fill_discrete(name = samVar) else
            scale_fill_gradient(name = samVar, low = "yellow", high = "blue")) + scale_shape_manual(values = shapeValues, name = shapeName)
}