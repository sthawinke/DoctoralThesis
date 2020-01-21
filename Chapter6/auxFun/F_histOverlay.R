#' Plot a histogram of z-values with normal overlay and vertical reference line
#' @param zValues the z-values to plot
#' @param breaks,xlab,main arguments passed on the hist() function
#' @param zSeq The support to plot
#' @param ... additional arguments passed on the hist() function
histOverlay = function(zValues, breaks = 30, xlab = "z-values",
                       zSeq = seq(min(zValues,-4), max(zValues, 4), by = 0.1),
                       main = "Histogram of z-values",
                       ylim = NULL,...){
  hist(zValues, breaks = breaks, xlab = xlab, freq = FALSE, ylim = ylim,
       main = main, xlim = range(zSeq),...)
  abline(v = 0, col = "red")
  lines(x = zSeq, y = dnorm(zSeq), col = "blue")
}