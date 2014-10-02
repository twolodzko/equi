\name{see}
\alias{see}
\alias{plot.see}
\title{Bootstrap standard error of equating (SEE)}
\description{
Compute bootstrap \emph{Standard Error of Equating} (SEE).
}
\usage{
see(data, fun, clusters, reps = 100)
\method{plot}{see}(x, type="l", col=1, lty=1, ylab="yx", xlab="x", alpha=0.5, \dots)
}
\arguments{
  \item{data}{a data frame.}
  \item{fun}{a function to pass to bootstrap. The function must return a concordance table with source test scores in first column named "\code{x}", and equated scores in the second column named "\code{yx}" as in \code{\link{equi}} or \code{\link{lin}} functions.}
  \item{reps}{number of bootstrap repetitions (100 by default).}
  \item{clusters}{vector or list of vectors of cluster groupping variables.}
  \item{x}{object to plot or print.}
  \item{type, col, lty, xlab, ylab}{arguments of the plot function.}
  \item{alpha}{plot lines opacity (range [0, 1]).}
  \item{\dots}{potentially further arguments passed from other methods.}
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
\code{\link{equi}}, \code{\link{smoothtab}}
}
\examples{
data(Tests)
x <- Tests[Tests$Sample == "P", "x"]
y <- Tests[Tests$Sample == "P", "y"]

data <- data.frame(x=x, y=y)

# a function to be passed to bootstrap

myfun1 <- function(data) {
	eq <- equi(smoothtab(data$x, data$y, presmoothing=TRUE))
	return(eq$conc)
}

# note: the number of iterations is small
#       only to run faster as an example

see(data, myfun1, reps=25)

## Bootstrap for NEAT-CE design

data(Tests)

myfun2 <- function(Tests) {
	p <- Tests[Tests$Sample == "P", 1:2]
	q <- Tests[Tests$Sample == "Q", 2:3]
	eq <- equi(smoothtab(p), smoothtab(q))
	return(eq$conc)
}

see(Tests, myfun2, reps=25)
}