\name{nslm}
\alias{nslm}
\alias{print.nslm}
\alias{summary.nslm}
\title{Function for predicting by natural cubic splines regression}
\description{
It uses Natural Cubic Splines linear regression for prediction.
}
\usage{
nslm(x, y, z = x, df = "BIC")
\method{print}{nslm}(x, \dots)
\method{summary}{nslm}(object, \dots)
}
\arguments{
  \item{x, y}{linear regression arguments, \code{y ~ a + bx}.}
  \item{z}{the values to be used in prediction from regression function.}
  \item{df}{the number of degrees of freedom for the Natural Cubic Spline used in \code{\link{nslm}} for estimating the missing values using AIC, BIC or the sum of squared errors. Possible values are: "BIC" (default), "AIC", "LSE" (for sum of squared errors) or numeric value for df given by the user.}
  \item{object}{an \code{summary} argument.}
  \item{\dots}{potentially further arguments passed from other methods.}
}

\value{
	\item{choice.method}{the method for choosing proper df, possible values: "none", "BIC", "AIC, "LSE".}
	\item{df}{number of degrees of freedom for \code{\link{ns}}.}
	\item{model}{linear regression model used.}
	\item{predicted}{predicted values.}
	\item{\dots}{potentially further arguments passed from other methods.}
}

\seealso{
\code{\link{ns}}, \code{\link{lm}}
}
\examples{
e <- runif(10, -1, 1)
plot(1:10, 1:10+e)
lines(nslm(1:10, 1:10+e, seq(1, 10, by=.2))$predicted~seq(1, 10, by=.2), type="l")

data(Tests)
x <- Tests[Tests$Sample == "P", "x"]
y <- Tests[Tests$Sample == "P", "y"]

(fit <- nslm(x, y))
summary(fit)
}

