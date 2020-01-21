#' A function to apply the Wilcoxon rank sum test or t-test for detecting differential absolute abundance
#' @param Y the matrix of sequencing counts
#' @param x The grouping factor
#' @param FC the total flow cytometry count vector
#' @param tieBreak a character vector, the tiebreak paradigm
#' @param S an estimate of the library sizes
#' @param test a function, which test to apply
wilcoxFC = function(tieBrokenY, x, test  = wilcox.test, statisticOnly = TRUE, nFac = NULL, seqX = NULL){
  apply(tieBrokenY, 2, function(z){if(statisticOnly) wilcox.test.fast(x = z[!x], y = z[x], nFac = nFac, seqX = seqX) else test(z ~ x)})
}
wilcox.test.fast = function(x, y, nFac, seqX){
  sum(rank(c(x, y))[seqX]) - nFac
}
