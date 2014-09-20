\name{trun}
\alias{trun}
\title{Truncate}
\description{
Truncate a vector given lower and upper margins.
}
\usage{
trun(x, lower, upper, margin=0)
}
\arguments{
  \item{x}{numeric vector.}
  \item{lower}{lower margin.}
  \item{upper}{upper margin.}
  \item{margin}{value of \code{margin} is subtracted form lower value of \code{range} and added to upper value of \code{range}.}
  }
\examples{
trun(1:10, 4, 8, margin=0)
trun(1:10, 4, 8, margin=0.5)
}
