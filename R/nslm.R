

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
	
	UseMethod("nslm")
}


print.nslm <- function(x, ...) {
	print(x$model, ...)
}


summary.nslm <- function(object, ...) {
	cat(paste("Natural Cubic Spline linear regression with",
						object$df, "degrees of freedom"))
	if (object$choice.method != "none")
		cat(paste(" choosen using", object$choice.method))
	cat(".\n")
	print(object$model, ...)
	print(summary(object$model), ...)
}

