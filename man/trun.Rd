\name{trun}
\alias{trun}
\title{Truncate}
\description{
Truncate a vector given a \code{range}.
}
\usage{
trun(x, range, margin = TRUE)
}
\arguments{
  \item{x}{numeric vector.}
  \item{range}{two-item vector containing lower and upper margin values.}
  \item{margin}{value of \code{margin} is subtracted form lower value of \code{range} and added to upper value of \code{range}.}
  }
\examples{
trun(1:10, c(4, 8), margin=TRUE)
trun(1:10, c(4, 8), margin=FALSE)
}
