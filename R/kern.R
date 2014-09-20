

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


print.kernsmooth <- function(x, ...) {
	print(x$data, ...)
	cat(paste("\nKernel Smoothing bandwidth:",
						round(x$h, digits=getOption("digits"))))
}


plot.kernsmooth <- function(x, type="l", add=FALSE, ...) {
	if (!add) {
		plot(x$data, type=type, ...)
	} else lines(x$data, type=type, ...)
}


cdfplot.kernsmooth <- function(x, type="s", add=FALSE, ...) {
	if (!add) {
		plot(x$data$x, cumsum(x$data$Freq)/sum(x$data$Freq), 
				 type=type, ylab="prob", xlab="", ...)
	} else lines(x$data$x, cumsum(x$data$Freq)/sum(x$data$Freq), 
							 type=type, ...)
}

