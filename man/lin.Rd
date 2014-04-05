\name{lin}
\alias{lin}
\title{
Linear equating
}
\description{
Function for linear equating
}
\usage{
lin(x, y, m=mean(y, na.rm=TRUE),
    s=sd(y, na.rm=TRUE))
}
\arguments{
  \item{x, y}{numerical vectors. Vector \code{y} can be omitted if \code{m} and \code{s} are provided.}
  \item{m}{target mean value, by default computed from \code{y}.}
  \item{s}{target sd value, by default computed from \code{y}.}
}
\value{
Returns object of class \code{lin}, that can be used for transforming scores using linear equating.
}
\references{
Kolen, M.J. & Brennan, R.J. (2004). \emph{Test Equating, Scaling, and Linking: Methods and Practices.} New York: Springer-Verlag.
}
\examples{
data(Tests)
x <- Tests[Tests$Sample == "P", "x"]
y <- Tests[Tests$Sample == "P", "y"]

(leq <- lin(x, y))
lin(x, leq)

lin(x, m=100, s=15)
}
