\name{smoothtab}
\alias{smoothtab}
\title{Presmoothing and postsmoothing of empirical distributions}
\description{
Function for \emph{Log-linear} presmoothing and/or \emph{Gaussian Kernel} postsmoothing.
}
\usage{
smoothtab(x, y, presmoothing = FALSE, postsmoothing = FALSE,
          bandwidth = 0.33, grid = 100, lldeg = 4,
          llxdeg = 1, raw = TRUE, cdf = TRUE)
}
\arguments{
  \item{x, y}{numeric vectors.}
  \item{presmoothing}{if \code{TRUE} \emph{Log-linear} presmoothing is applied.}
  \item{postsmoothing}{if \code{TRUE} \emph{Gaussian Kernel} postsmoothing is applied.}
  \item{bandwidth}{sets bandwidth for \emph{Kernel Smoothing}. Default value is \code{0.33}. If set to \code{"auto"}, it uses \code{\link[KernSmooth]{dpik}} for automatic selection of bandwidth.
  }
  \item{grid}{sets \code{gridsize} for \code{\link{bkde}} or \code{\link{bkde2D}}.}
  \item{lldeg}{degree of the polynomial in log-linear presmoothing (deafult is \code{4}).}
  \item{llxdeg}{degree of the polynomial in log-linear presmoothing for interaction term (deafult is \code{1}).}
  \item{raw}{if \code{TRUE} computes raw polynomials for log-linear presmoothing, see: \code{\link[splines]{poly}}}
  \item{cdf}{if \code{FALSE} compute probabilities rather than cumulative probabilities (default).}
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
\note{See also \code{\pkg{equate}} and \code{\pkg{kequate}} packages.}
\seealso{
\code{\link{equi}}, \code{\link[splines]{poly}}, \code{\link[KernSmooth]{bkde}}, \code{\link[KernSmooth]{bkde2D}}
}
\examples{
data(Tests)
x <- Tests[Tests$Sample == "P", "x"]
y <- Tests[Tests$Sample == "P", "y"]
smoothtab(x)
smoothtab(x, presmoothing=TRUE)
smoothtab(x, y, presmoothing=TRUE, postsmoothing=TRUE)
}