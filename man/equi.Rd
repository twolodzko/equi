% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/equi.R
\name{equi}
\alias{equi}
\alias{plot.equi}
\title{Equipercentile equating}
\usage{
equi(x, y, range, truncate = TRUE, df = "BIC", margin = 0.5)

\method{plot}{equi}(x, diff = FALSE, ref = TRUE, type = "l", ...)
}
\arguments{
\item{x}{smoothtab object or numeric vector. For NEAT first column of smoothtab
data frame is source test and second is anchor test (see Details).}

\item{y}{smoothtab or equi object. For NEAT first column of smoothtab data frame
is anchor test and second is target test (see Details).}

\item{range}{two-item vector of minimum and maximum score points of \code{y}.}

\item{truncate}{if \code{TRUE} truncates the values using trun.}

\item{df}{degrees of freedom for \code{ns} splines. If \code{df} is set to "BIC", "AIC" or "LSE"
it uses BIC, AIC or sum of squared errors for choosing the best df value.}

\item{margin}{see: \code{\link{trun}}.}

\item{diff}{plot a difference from identity function.}

\item{ref}{plot a reference line, i.e. an identity function.}

\item{type}{type of the plot.}

\item{\dots}{potentially further arguments passed from other methods.}
}
\value{
Returns object of class equi, that can be used for transforming test scores using
equi function (see examples).
}
\description{
Functions for equipercentile equating in Equivalent Groups (EG), Single Group (SG)
and Nonequivalent groups with Anchor Test using Chained Equating (NEAT-CE) designs.
}
\details{
If \code{x} is smoothtab object and \code{y} is not provided, single group equating is applied.
If \code{x} and \code{y} are smoothtab objects, equivalent groups or NEAT-CE equating is applied,
depending on the kind of data provided. If \code{x} is numeric vector and \code{y} is equi object,
\code{x} is transformed using equating function provided in \code{y}.
}
\examples{
# Single Group design with presmoothing

data(Tests)
x <- Tests[Tests$Sample == "P", "x"]
y <- Tests[Tests$Sample == "P", "y"]

(eq <- equi(smoothtab(x, y, presmoothing=TRUE)))
yx <- equi(x, eq)

moments(list(x, y, yx))

# Equivalent Groups design

data(Tests)
x <- Tests[Tests$Sample == "P", "x"]
y <- Tests[Tests$Sample == "P", "y"]

(eq <- equi(smoothtab(x), smoothtab(y)))
yx <- equi(x, eq)

# NEAT-CE design, 'x' is equated to 'y' via anchor test 'a' (xa, ya).
# in 'p' first column is source test, second is anchor test
# in 'q' first column is anchor test, second is target test

data(Tests)
p <- Tests[Tests$Sample == "P", 1:2]
q <- Tests[Tests$Sample == "Q", 2:3]

(eq <- equi(smoothtab(p), smoothtab(q)))
yx <- equi(p[, 1], eq)
}
\references{
Kolen, M.J. & Brennan, R.J. (2004). \emph{Test Equating, Scaling, and Linking:
Methods and Practices.} New York: Springer-Verlag.

von Davier, A.A., Holland, P.W. & Thayer, D.T. (2004). \emph{The Kernel Method
of Test Equating.} New York: Springer-Verlag.

Green, P.J. & Silverman, B.W. (1993). \emph{Nonparametric Regression and Generalized
Linear Models: A roughness penalty approach.} London: Chapman & Hall/CRC.
}
\seealso{
\code{\link{ns}}, \code{\link{bs}}, \code{\link{poly}}, \code{\link{approx}}
}

