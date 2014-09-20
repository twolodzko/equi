\name{cdfplot}
\alias{cdfplot}
\alias{cdfplot.default}
\title{Plot empirical cumulative distribution function}
\description{
Creates empirical cumulative distribution function plot.
}
\usage{
\method{cdfplot}{default}(x, add=FALSE, ylab="F(x)", xlab="",
        lwd=1, main="", \dots)
}
\arguments{
  \item{x}{numeric vector.}
  \item{add}{if \code{TRUE} adds new plot to existing one.}
  \item{ylab, xlab, lwd, main, \dots}{pass additional arguments to \code{\link[graphics]{plot}}. See also: \code{\link[graphics]{plot.default}}}
  }
\value{
Returns a table with empirical cumulative probabilities.
}
\references{
Cai, E. (June 25, 2013). \emph{Exploratory Data Analysis: 2 Ways of Plotting Empirical Cumulative Distribution Functions in R.} \href{http://chemicalstatistician.wordpress.com/2013/06/25/exploratory-data-analysis-2-ways-of-plotting-empirical-cumulative-distribution-functions-in-r/}{The Chemical Statistician blog.}
}
\examples{
cdfplot(rpois(1000, 5))
}