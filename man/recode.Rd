\name{recode}
\alias{recode}
\title{
Recode variable
}
\description{
Recode a vector \code{from} chosen values \code{to} new values. Any value not existing in \code{from} vector will be recoded to \code{other}. Numeric vector can be truncated.
}
\usage{
recode(x, from, to, other, truncate = TRUE, margin=0.5)
}
\arguments{
  \item{x}{vector.}
  \item{from}{source values.}
  \item{to}{target values.}
  \item{other}{if provided, other values are recoded to \code{other}.}
  \item{truncate}{if \code{TRUE} values below \code{\link{min}} and above \code{\link{max}} of the \code{to} vector are recoded to \code{\link{min}} and \code{\link{max}} respectively. See: \code{\link{trun}}}
  \item{margin}{see: \code{\link{trun}}}
}
\seealso{
\code{\link{trun}}
}
\examples{
recode(1:10, from=3:5, to=5:3, truncate=TRUE)
recode(1:10, from=3:5, to=5:3, other=NA)
}