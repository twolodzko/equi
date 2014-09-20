\name{kernsmooth}
\alias{kernsmooth}
\alias{print.kernsmooth}
\alias{plot.kernsmooth}
\alias{cdfplot.kernsmooth}
\title{
Kernel Smoothing
}
\description{
Kernel Smoothing as described by von Davier, Holland & Thayer (2004). It can be used for demonstrational pourposes as lower level function is used internally in \code{\link{smoothtab}} for postsmoothing and there is no need for using \code{kernsmooth} directly.
}
\usage{
kernsmooth(x, h = "auto", grid = 100)
\method{print}{kernsmooth}(x, \dots)
\method{plot}{kernsmooth}(x, type="l", add=FALSE, \dots)
\method{cdfplot}{kernsmooth}(x, type="s", add=FALSE, \dots)
}

\arguments{
  \item{x}{vector of values.}
  \item{h}{smoothing bandwidth: as h goes bigger the distribution is transformed to into normal distribution. If set to "auto" (default) it uses automatic algorithm for selecting its value (von Davier, Holland & Thayer, 2004).}
  \item{grid}{number of points to be created.}
  \item{type}{type of the plot.}
  \item{add}{add a plot to previous one.}
  \item{\dots}{potentially further arguments passed from other methods.}
}

\references{
von Davier, A.A., Holland, P.W. & Thayer, D.T. (2004). \emph{The Kernel Method of Test Equating.} New York: Springer-Verlag.
}

\examples{
data(Tests)
x <- Tests[Tests$Sample == "P", "x"]

ct <- conttab(x, prob=TRUE)

plot(ct, col="red", ylim=c(0, 0.13))

for (i in seq(.35, 2.5, by=.1))
  lines(kernsmooth(x, h=i, grid=1000)$data)

(ks <- kernsmooth(x))
plot(ks)
cdfplot(ks)
}

