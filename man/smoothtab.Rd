\name{smoothtab}
\alias{smoothtab}
\alias{is.smoothtab}
\alias{as.smoothtab}
\alias{print.smoothtab}
\alias{plot.smoothtab}
\alias{cdfplot.smoothtab}
\title{Presmoothing and postsmoothing of empirical distributions}
\description{
Function for \emph{Log-linear} presmoothing and/or \emph{Gaussian Kernel} postsmoothing.
}
\usage{
smoothtab(x, y, presmoothing=FALSE, postsmoothing=FALSE,
          bandwidth="auto", lldeg=4, llxdeg=1,
          raw=TRUE, cdf=TRUE, margin=0.5, grid=100)
is.smoothtab(x)
as.smoothtab(x)
\method{print}{smoothtab}(x, \dots)
\method{plot}{smoothtab}(x, type="s", lty=1:6, add=FALSE, \dots)
\method{cdfplot}{smoothtab}(x, add=FALSE, \dots)
}
\arguments{
  \item{x, y}{numeric vectors.}
  \item{presmoothing}{if \code{TRUE} \emph{Log-linear} presmoothing is applied.}
  \item{postsmoothing}{if \code{TRUE} \emph{Gaussian Kernel} postsmoothing is applied.}
  \item{bandwidth}{sets bandwidth for \emph{Kernel Smoothing}. Use \code{"auto"} (default) for automatic selection of bandwidth.
  }
  \item{lldeg}{degree of the polynomial in log-linear presmoothing (deafult is \code{4}).}
  \item{llxdeg}{degree of the polynomial in log-linear presmoothing for interaction term (deafult is \code{1}).}
  \item{raw}{if \code{TRUE} computes raw polynomials for log-linear presmoothing, see: \code{\link[stats]{poly}}}
  \item{cdf}{if \code{FALSE} compute probabilities rather than cumulative probabilities (default).}
  \item{margin}{if \code{postsmoothing=TRUE}, it defines the margins of points range to be created using \emph{Gaussian Kernel} postsmoothing. The function returns point values +/- margin.}
  \item{grid}{if \code{postsmoothing=TRUE}, it defines the number of points to be created using \emph{Gaussian Kernel} postsmoothing.}
  \item{type}{type of the plot.}
  \item{add}{add a plot to previous one.}
  \item{lty}{a vector of line types, see \code{\link{par}}.}
  \item{\dots}{potentially further arguments passed from other methods.}
}
\value{
Returns two-column \code{data.frame} with unique score points and coresponding probabilities, or a list of two \code{data.frames} for marginal probabilities of joint distributions.
}
\references{
Holland, P.W. & Thayer, D.T. (2000). Univariate and Bivariate Loglinear Models for Discrete Test
Score Distributions. \emph{Journal of Educational and Behavioral Statistics, 25}(2), 133-183.

Kolen, M.J. & Brennan, R.J. (2004). \emph{Test Equating, Scaling, and Linking: Methods and Practices.} New York: Springer-Verlag.

von Davier, A.A., Holland, P.W. & Thayer, D.T. (2004). \emph{The Kernel Method of Test Equating.} New York: Springer-Verlag.

Wand, M.P. & Jones, M.C. (1995). \emph{Kernel Smoothing.} London: Chapman & Hall/CRC. 
}
\note{See also \pkg{equate} and \pkg{kequate} packages.}
\examples{
data(Tests)
x <- Tests[Tests$Sample == "P", "x"]
y <- Tests[Tests$Sample == "P", "y"]

(st <- smoothtab(x))
plot(st)

smoothtab(x, presmoothing=TRUE)
smoothtab(x, y, presmoothing=TRUE, postsmoothing=TRUE)
}