

#' Function for predicting by natural cubic splines regression
#' 
#' It uses Natural Cubic Splines linear regression for prediction. 
#' 
#' @param x,y linear regression arguments, \code{y ~ a + bx}.
#' @param z   the values to be used in prediction from regression function.
#' @param df  the number of degrees of freedom for the Natural Cubic Spline
#'            used in \code{nslm} for estimating the missing values using AIC, BIC
#'            or the sum of squared errors. Possible values are: "BIC" (default),
#'            "AIC", "LSE" (for sum of squared errors) or numeric value for df
#'            given by the user.
#' 
#' @return
#' 
#' choice.method the method for choosing proper df, possible values:
#'               "none", "BIC", "AIC, "LSE".
#' df            number of degrees of freedom for ns.
#' model         linear regression model used.
#' predicted 	   predicted values.
#' 
#' @seealso \code{\link{ns}}, \code{\link{lm}}
#' 
#' @examples
#' 
#' e <- runif(10, -1, 1)
#' plot(1:10, 1:10+e)
#' lines(nslm(1:10, 1:10+e, seq(1, 10, by=.2))$predicted~seq(1, 10, by=.2), type="l")
#' 
#' data(Tests)
#' x <- Tests[Tests$Sample == "P", "x"]
#' y <- Tests[Tests$Sample == "P", "y"]
#' 
#' (fit <- nslm(x, y))
#' summary(fit)
#' 
#' @export

nslm <- function(x, y, z=x, df="BIC") {
	out <- list()
	if (df %in% c("AIC", "BIC", "LSE")) {
		k <- length(unique(y))
		crit <- NULL
		for (i in 1:(k-1)) {
			fit <- lm(y~ns(x, df=i))
			if (df == "LSE") {
				pred <- predict(fit, data.frame(x=x))
				crit[i] <- sum((y-pred)^2)
			} else if (df == "AIC") { crit[i] <- AIC(fit)
			} else crit[i] <- BIC(fit)
		}
		out$choice.method <- df
		out$criteria <- crit
		out$df <- which.min(crit)
	} else {
		out$choice.method <- "none"
		out$df <- df
	}
	out$model <- lm(y~ns(x, df=out$df))
	out$predicted <- predict(out$model, data.frame(x=z))
	class(out) <- "nslm"
	
	return(out)
}


#' @export

print.nslm <- function(x, ...) {
	print(x$model, ...)
}


#' @rdname nslm
#' 
#' @param object  an summary argument.
#' @param \dots  potentially further arguments passed from other methods.
#' 
#' @export

summary.nslm <- function(object, ...) {
	cat(paste("Natural Cubic Spline linear regression with",
						object$df, "degrees of freedom"))
	if (object$choice.method != "none")
		cat(paste(" choosen using", object$choice.method))
	cat(".\n")
	print(object$model, ...)
	print(summary(object$model), ...)
}

