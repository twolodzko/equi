

#' Plot empirical cumulative distribution function
#' 
#' Creates empirical cumulative distribution function plot.
#' 
#' @param x                        numeric vector.
#' @param add                      if \code{TRUE} adds new plot to existing one.
#' @param ylab,xlab,lwd,main,\dots pass additional arguments to \code{\link[graphics]{plot}}.
#'                                 See also: \code{\link{plot.default}}.
#' 
#' @return
#' Returns a table with empirical cumulative probabilities.
#' 
#' @references
#' Cai, E. (June 25, 2013). \emph{Exploratory Data Analysis: 2 Ways of Plotting Empirical Cumulative Distribution Functions in R.}
#' \href{http://chemicalstatistician.wordpress.com/2013/06/25/exploratory-data-analysis-2-ways-of-plotting-empirical-cumulative-distribution-functions-in-r/}{The Chemical Statistician blog.}
#' 
#' @examples
#' 
#' cdfplot(rpois(1000, 5))
#' 
#' @export

cdfplot <- function(x, ...) UseMethod("cdfplot")


#' @export

cdfplot.default <- function(x, add=FALSE, ylab="F(x)", xlab="",
                            lwd=1, main="", ...) {
  x <- sort(x)
  n <- length(x)
  nx <- (1:n)/n
  if (!add) {
    plot(x, nx, ylab=ylab, xlab=xlab,
         type="s", lwd=lwd, main=main, ...)
  } else lines(x, nx, lwd=lwd, ...)
  invisible(tapply(nx, x, sum))
}

