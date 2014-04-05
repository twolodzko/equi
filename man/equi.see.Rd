\name{equi.see}
\alias{equi.see}
\title{Bootstrap Standard Error of Equating (SEE)}
\description{
Compute bootstrap \emph{Standard Error of Equating} (SEE).
}
\usage{
equi.see(x, y, design = "sg", reps = 100,
        presmoothing = FALSE, postsmoothing = FALSE,
        fun = "ns", df = 10, bandwidth = 0.33,
        truncate = TRUE, range, range.y, grid = 100,
        lldeg = 4, llxdeg = 1, raw = TRUE, clusters)
}
\arguments{
  \item{x, y}{numeric vectors or two-column marices.}
  \item{design}{desing function: \code{"sg"} for Single Group, \code{"eg"} for Equivalent Groups (default) and \code{"neat"} for Nonequivalent groups with Anchor Test.}
  \item{reps}{number of bootstrap repetitions (\code{100} by default).}
  \item{clusters}{vector or list of vectors of cluster groupping variables.}
  \item{range, range.y}{ranges for truncating x (and y for NEAT). See: \code{\link{equi}}, \code{\link{trun}}. }
  \item{presmoothing, postsmoothing, fun, df, truncate}{see: \code{\link{equi}}}
  \item{bandwidth, grid, lldeg, llxdeg, raw}{see: \code{\link{smoothtab}}}
}
\value{
Returns \emph{Bias}, \emph{SE} and \emph{RMSE} values for each score point and averaged SEE values coputed with mean or average weighted on score points probabilities.
}
\references{
Efron, B. & Tibshirani, R.J. (1993). \emph{An Introduction to the Bootstrap.} London: Chapman & Hall/CRC.

Field, C.A. & Welsh, A.H. (2007). Bootstrapping clustered data. \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology), 69}(3), 369-390.

Kolen, M.J. & Brennan, R.J. (2004). \emph{Test Equating, Scaling, and Linking: Methods and Practices.} New York: Springer-Verlag.

Rena, S., Lai, H., Tong, W., Aminzadeh, M., Hou, X. & Lai, S. (2010). Nonparametric bootstrapping for hierarchical data. \emph{Journal of Applied Statistics, 37}(9), 1487-1498.

von Davier, A.A., Holland, P.W. & Thayer, D.T. (2004). \emph{The Kernel Method of Test Equating.} New York: Springer-Verlag.

Wang, C. (2011). \emph{An investigation of bootstrap methods for estimating the standard error of equating under the common-item nonequivalent groups design.} doctoral PhD diss., University of Iowa. \url{http://ir.uiowa.edu/etd/1188}
}
\seealso{
\code{\link{equi}}, \code{\link{equi.ce}}, \code{\link{smoothtab}}
}
\examples{
data(Tests)
x <- Tests[Tests$Sample == "P", "x"]
y <- Tests[Tests$Sample == "P", "y"]

equi.see(x, y, design="sg", reps=10)
equi.see(x, y, design="eg", reps=10)

data(Tests)
p <- Tests[Tests$Sample == "P", 1:2]
q <- Tests[Tests$Sample == "Q", 2:3]

equi.see(p, q, design="neat", reps=10)

# SEE for clustered data

data(Tests)
x <- Tests[Tests$Sample == "P", "x"]
y <- Tests[Tests$Sample == "P", "y"]
z <- Tests[Tests$Sample == "P", "School"]

equi.see(x, y, design="sg", clusters=z, reps=10)
equi.see(x, y, design="eg", clusters=list(z, z), reps=10)

data(Tests)
p <- Tests[Tests$Sample == "P", c(1:2, 5)]
q <- Tests[Tests$Sample == "Q", c(2:3, 5)]

equi.see(p, q, design="neat", clusters=list(p[, 3], q[,3]), reps=10)
}