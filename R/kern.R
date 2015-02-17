

kern <- function(x, xj, p, m, s, h="auto", hmin=0.1, hmax=1, k=1, pen=FALSE) {
	if (h == "auto") {
		f <- function(h) kern(x, xj, p, m, s, h, k=k, pen=TRUE)$pen
		return(kern(x, xj, p, m, s, h=optimize(f, lower=hmin, upper=hmax)$minimum))
	} else {
		k <- length(xj)
		n <- length(x)
		a <- sqrt(s/(s+h^2))
		res <- NULL
		
		Rjx <- function(x, a, xj, m, h, i)
			(x[i] - a*xj - (1-a)*m) / (a*h)
		
		for (i in 1:n) {
			res[i] <- sum(p * dnorm(Rjx(x, a, xj, m, h, i)) / (a*h))
		}
		
		out <- list(data=data.frame(x=x, Freq=res), h=h)
		class(out) <- "kernsmooth"
		
		if (pen) {
			p.hat <- NULL
			der <- NULL
			for (i in 1:k) {
				p.hat[i] <- sum(p * dnorm(Rjx(xj, a, xj, m, h, i)) / (a*h))
				der[i] <- -sum(p * dnorm(Rjx(x, a, xj, m, h, i)) / (a*h) )
			}
			p1 <- sum((p-p.hat)^2)
			p2 <- sum( (der > 0) * (1 - (der <= 0) ))
			out$pen <- p1+k*p2
		} 
		return(out)
	}
}


#' Kernel smoothing 
#' 
#' Kernel Smoothing as described by von Davier, Holland & Thayer (2004). It can be used
#' for demonstrational pourposes while a lower level function is used internally in smoothtab
#' for postsmoothing and there is no need for using kernsmooth directly. 
#' 
#' @param x    vector of values.
#' @param h    smoothing bandwidth: as h goes bigger the distribution is transformed to into
#'             normal distribution. If set to "auto" (default) it uses automatic algorithm for
#'             selecting its value (von Davier, Holland & Thayer, 2004).
#' @param grid number of points to be created.
#' 
#' @references
#' 
#' von Davier, A.A., Holland, P.W. & Thayer, D.T. (2004). \emph{The Kernel Method
#' of Test Equating.} New York: Springer-Verlag.
#' 
#' @examples
#' 
#' data(Tests)
#' x <- Tests[Tests$Sample == "P", "x"]
#' 
#' ct <- conttab(x, prob=TRUE)
#' 
#' plot(ct, col="red", ylim=c(0, 0.13))
#' 
#' for (i in seq(.35, 2.5, by=.1))
#'   lines(kernsmooth(x, h=i, grid=1000)$data)
#' 
#' (ks <- kernsmooth(x))
#' plot(ks)
#' cdfplot(ks)
#' 
#' @export

kernsmooth <- function(x, h="auto", grid=100) {
	x <- x[!is.na(x)]
	m <- mean(x)
	s <- sd(x)
	ct <- conttab(x)
	xj <- ct$x
	p <- ct$Freq/sum(ct$Freq)
	x.new <- seq(min(xj), max(xj), length.out=grid)
	return(kern(x.new, xj, p, m, s, h=h))
}


#' @export

print.kernsmooth <- function(x, ...) {
	print(x$data, ...)
	cat(paste("\nKernel Smoothing bandwidth:",
            round(x$h, digits=getOption("digits"))))
}


#' @rdname kernsmooth
#' 
#' @param type  type of the plot.
#' @param add   add a plot to previous one.
#' @param \dots potentially further arguments passed from other methods.
#' 
#' @export

plot.kernsmooth <- function(x, type="l", add=FALSE, ...) {
	if (!add) {
		plot(x$data, type=type, ...)
	} else lines(x$data, type=type, ...)
}


#' @rdname kernsmooth
#' @export

cdfplot.kernsmooth <- function(x, type="s", add=FALSE, ...) {
	if (!add) {
		plot(x$data$x, cumsum(x$data$Freq)/sum(x$data$Freq),
         type=type, ylab="prob", xlab="", ...)
	} else lines(x$data$x, cumsum(x$data$Freq)/sum(x$data$Freq),
               type=type, ...)
}

