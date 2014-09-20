

see <- function(data, fun, clusters, reps = 100) {
	
	n <- nrow(data)
	tmp <- do.call(fun, list(data))
	yx <- tmp$yx
	k <- length(tmp$x)
	cmp <- matrix(NA, reps, k)
	
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
		rownames(tmp) <- tmp$x
		cmp[i, ] <- tmp[as.character(tmp$x), "yx"]
		
		setTxtProgressBar(pb, i)
	}
	close(pb)
	
	bias <- colMeans(cmp, na.rm = TRUE) - yx
	se <- apply(cmp, 2, sd, na.rm = TRUE)
	rmse <- sqrt(bias^2 + se^2)
	
	out <- data.frame(bias=bias, se=se, rmse=rmse)
	
	print(colMeans(out))
	invisible(out)
}

