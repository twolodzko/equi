\name{moments}
\alias{moments}
\title{Moments}
\description{Calculates moments of given distribution up to the nth moment.}
\usage{
moments(x, num = 4, excess = TRUE)
}
\arguments{
  \item{x}{numeric vector or list of vectors.}
  \item{num}{calculate up to the nth moment given by \code{num}.}
  \item{excess}{if \code{TRUE} calculates excess kurtosis.}
}
\examples{
moments(rpois(1000, 5))
moments(rnorm(1000), num=6)
}