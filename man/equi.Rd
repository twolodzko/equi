\name{equi}
\alias{equi}
\title{
Equipercentile Equating
}
\description{
Functions for equipercentile equating in \emph{Equivalent Groups} (EG), \emph{Single Group} (SG) and \emph{Nonequivalent groups with Anchor Test} using \emph{Chained Equating} (NEAT-CE) designs.
}
\usage{
equi(x, y, range, fun = "ns", df = "BIC",
     truncate = TRUE, margin=0.5)
}
\arguments{
  \item{x}{\code{\link{smoothtab}} object or numeric vector. For NEAT first column of \code{smoothtab} is source test and second is anchor test.}
  \item{y}{\code{\link{smoothtab}} or \code{\link{equi}} object. For NEAT first column of \code{smoothtab} is anchor test and second is target test.}
  \item{range}{two-item vector of minimum and maximum score points of \code{y} scale.}
  \item{xlim, ylim}{two-item vectors of minimum and maximum score points of \code{x} scale anchor test and \code{y} scale target test.}
  \item{fun}{functions used for interpolating equating result scores from equating functions. \code{"none"} uses \code{approx} function for linear interpolation, \code{"ns"} (default) and \code{"bs"} use linear regression with \emph{Natural Cubic Spline} and \emph{B-Spline} functions respectively. See: \code{\link[stats]{approx}}, \code{\link[splines]{ns}}, \code{\link[splines]{bs}}}
  \item{df}{degrees of freedom for \code{\link[splines]{ns}} splines. If \code{df} is set to "BIC", "AIC" or "LSE" it uses BIC, AIC or sum of squared errors for choosing the best \code{df} value. }
  \item{truncate}{if \code{TRUE} truncates the values using \code{\link{trun}.}}
  \item{margin}{see: \code{\link{trun}}}
}
\details{
If \code{x} is \code{smoothtab} object and \code{y} is not provided, single group equating is applied. If \code{x} and \code{y} are \code{smoothtab} objects, equivalent groups or NEAT-CE equating is applied, depending on the kind of data provided. If \code{x} is numeric vector and \code{y} is \code{equi} object, \code{x} is transformed using equating function provided in \code{y}.
}
\value{
Returns object of class \code{equi}, that can be used for transforming test scores using \code{equi} function (see examples).
}
\references{
Kolen, M.J. & Brennan, R.J. (2004). \emph{Test Equating, Scaling, and Linking: Methods and Practices.} New York: Springer-Verlag.

von Davier, A.A., Holland, P.W. & Thayer, D.T. (2004). \emph{The Kernel Method of Test Equating.} New York: Springer-Verlag.

Green, P.J. & Silverman, B.W. (1993). \emph{Nonparametric Regression and Generalized Linear Models: A roughness penalty approach.} London: Chapman & Hall/CRC.
}
\seealso{
\code{\link[splines]{ns}}, \code{\link[splines]{bs}}, \code{\link[splines]{poly}}, \code{\link[stats]{approx}}
}
\examples{
# Single Group design with presmoothing

data(Tests)
x <- Tests[Tests$Sample == "P", "x"]
y <- Tests[Tests$Sample == "P", "y"]

eq <- equi(smoothtab(x, y, presmoothing=TRUE))
yx <- equi(x, eq)

moments(list(x, y, yx))

# Equivalent Groups design

data(Tests)
x <- Tests[Tests$Sample == "P", "x"]
y <- Tests[Tests$Sample == "P", "y"]

eq <- equi(smoothtab(x), smoothtab(y))
yx <- equi(x, eq)

# NEAT-CE design, 'x' is equated to 'y' via anchor test 'a' (xa, ya).
# in 'p' first column is source test, second is anchor test
# in 'q' first column is anchor test, second is target test

data(Tests)
p <- Tests[Tests$Sample == "P", 1:2]
q <- Tests[Tests$Sample == "Q", 2:3]

eq <- equi(smoothtab(p), smoothtab(q))
yx <- equi(p[, 1], eq)
}
