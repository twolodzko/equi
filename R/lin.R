

#' Linear equating 
#' 
#' Function for linear equating 
#' 
#' @param x,y numerical vectors. Vector \code{y} can be omitted if
#'            \code{m} and \code{s} are provided.
#' @param m   target mean value, by default computed from \code{y}.
#' @param s   target sd value, by default computed from \code{y}.
#' 
#' @return
#' 
#' Returns object of class lin, that can be used for transforming
#' scores using linear equating. 
#' 
#' @references
#' 
#' Kolen, M.J. & Brennan, R.J. (2004). \emph{Test Equating, Scaling, and Linking:
#' Methods and Practices.} New York: Springer-Verlag. 
#' 
#' @examples
#' 
#' data(Tests)
#' x <- Tests[Tests$Sample == "P", "x"]
#' y <- Tests[Tests$Sample == "P", "y"]
#' 
#' (leq <- lin(x, y))
#' summary(leq)
#' lin(x, leq)
#' 
#' lin(x, m=100, s=15)
#' 
#' @export

lin <- function(x, y, m=mean(y, na.rm=TRUE), s=sd(y, na.rm=TRUE)) {
	
	if (!missing(y) && class(y) == "lin")
		return(x * y$slope + y$intercept)
	
	slope <- s/sd(x, na.rm=TRUE)
	intercept <- m - slope * mean(x, na.rm=TRUE)
	
	out <- list(slope=slope, intercept=intercept)
	class(out) <- c("lin")
	out$yx <- lin(x, out)
	
	return(out)
}


#' @export

print.lin <- function(x, ...){
	return(print(data.frame(slope=x$slope, intercept=x$intercept), ...))
}


#' @rdname lin
#' @export

summary.lin <- function(object, ...) {
	cat(paste("y = ", round(object$intercept, digits=getOption("digits")), " + ",
            round(object$slope, digits=getOption("digits")), " * x", sep=""))
}

