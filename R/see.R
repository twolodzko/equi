

see <- function(data, fun, clusters, reps = 100) {
	
	n <- nrow(data)
	est <- do.call(fun, list(data))
	k <- length(est$x)
	cmp <- matrix(NA, reps, k)
	colnames(cmp) <- est$x
	
	pb <- txtProgressBar(min=1, max=reps, initial=0, width=50, style=3)
	
	for (i in 1:reps) {
		if (missing(clusters)) {
			index <- sample(1:n, n, replace = TRUE)
		}
		else {
			unq <- unique(clusters)
			nu <- length(unq)
			clsamp <- sample(unq, nu, replace = TRUE)
			index <- NULL
			for (j in clsamp) {
				index <- c(index, which(clusters == j))
			}
		}
		samp <- data[index, ]
		
		tmp <- do.call(fun, list(samp))
		tmp.yx <- tmp$yx
		names(tmp.yx) <- tmp$x
		cmp[i, ] <- tmp.yx[as.character(est$x)]
		
		setTxtProgressBar(pb, i)
	}
	close(pb)
	
	bias <- colMeans(cmp, na.rm = TRUE) - est$yx
	se <- apply(cmp, 2, sd, na.rm = TRUE)
	rmse <- sqrt(bias^2 + se^2)
	
	out <- list(bootstrapSample=cmp,
							estimate=est,
							see=data.frame(bias=bias, se=se, rmse=rmse))
	class(out) <- "see"
	return(out)
}


print.see <- function(x, ...) {
	print(colMeans(x$see))
}


plot.see <- function(x, type="l", col=1, lty=1, ylab="yx", xlab="x", alpha=0.5, ...) {
	col <- col2rgb(col)/255
	col <- rgb(col[1], col[2], col[3], alpha=alpha)
	matplot(t(x$bootstrapSample), type=type, col=col, lty=lty, ylab=ylab, xlab=xlab, ...)
	lines(x$estimate$yx, col="red", lwd=3)
}

