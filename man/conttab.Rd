\name{conttab}
\alias{conttab}
\title{
Create contingency table
}
\usage{
conttab(..., numeric = TRUE)
}
\arguments{
  \item{\dots}{contingency table is build of variables in \code{\dots}}
  \item{numeric}{if \code{TRUE} (default) variables in contingency table are converted to \code{numeric} variables, else thay are treated as \code{factor}.}
}

\seealso{
\code{\link{table}}
}
\examples{
data(Tests)
x <- Tests[Tests$Sample == "P", "x"]
y <- Tests[Tests$Sample == "P", "y"]

conttab(x)
conttab(x, y)
}