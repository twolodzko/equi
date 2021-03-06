% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/lin.R
\name{lin}
\alias{lin}
\alias{summary.lin}
\title{Linear equating}
\usage{
lin(x, y, m = mean(y, na.rm = TRUE), s = sd(y, na.rm = TRUE))

\method{summary}{lin}(object, ...)
}
\arguments{
\item{x,y}{numerical vectors. Vector \code{y} can be omitted if
\code{m} and \code{s} are provided.}

\item{m}{target mean value, by default computed from \code{y}.}

\item{s}{target sd value, by default computed from \code{y}.}
}
\value{
Returns object of class lin, that can be used for transforming
scores using linear equating.
}
\description{
Function for linear equating
}
\examples{
data(Tests)
x <- Tests[Tests$Sample == "P", "x"]
y <- Tests[Tests$Sample == "P", "y"]

(leq <- lin(x, y))
summary(leq)
lin(x, leq)

lin(x, m=100, s=15)
}
\references{
Kolen, M.J. & Brennan, R.J. (2004). \emph{Test Equating, Scaling, and Linking:
Methods and Practices.} New York: Springer-Verlag.
}

