

lin <- function(x, y, m=mean(y, na.rm=TRUE),
												s=sd(y, na.rm=TRUE)) {
	
	if (!missing(y) && class(y) == "lin")
		return(x * y$slope + y$intercept)
	
	slope <- s/sd(x, na.rm=TRUE)
	intercept <- m - slope * mean(x, na.rm=TRUE)
	
	out <- list(slope=slope, intercept=intercept)
	class(out) <- c("lin")
	out$yx <- lin(x, out)
	
	return(out)
}


print.lin <- function(x, ...){
	return(print(data.frame(slope=x$slope, intercept=x$intercept), ...))
}


summary.lin <- function(object, ...) {
	cat(paste("y = ", round(object$intercept, digits=getOption("digits")), " + ",
						round(object$slope, digits=getOption("digits")), " * x", sep=""))
}

