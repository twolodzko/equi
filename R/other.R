

recode <- function(x, from, to, other, truncate=TRUE, margin=0.5) {
  new_x <- x
  k <- length(from)
  range <- c(min(to, na.rm=TRUE), max(to, na.rm=TRUE))
  for (i in 1:k) new_x[x == from[i]] <- to[i]
  if (!missing(other)) new_x[!(x %in% from)] <- other
  if (truncate) new_x <- trun(new_x, lower=range[1], upper=range[2], margin=margin)
  return(new_x)
}


trun <- function(x, lower, upper, margin=0) {
  if (!missing(lower)) {
  	lower <- lower - margin
  	x[x < lower] <- lower
  }
  if (!missing(upper)) {
  	upper <- upper + margin
  	x[x > upper] <- upper
  }
  return(x)
}


moments <- function(x, num=4, excess=TRUE) {
  if (is.list(x))
    return(sapply(x, moments, num=num, excess=excess))
  x <- x[!is.na(x)]
  m <- mean(x)
  s <- sd(x)
  out <- array(c(m, s))
  if (num > 2) out[3] <- mean(((x-m)/s)^3)
  if (num > 3) {
    out[4] <- mean(((x-m)/s)^4)
    if (excess) out[4] <- out[4] - 3
  }
  if (num > 4) for (i in 5:num) out[i] <- mean((x-m)^i)
  names(out) <- c("mean", "sd", "skew", "kurt", 5:num)[1:num] 
  return(out)
}


conttab <- function(..., numeric=TRUE, prob=FALSE) {
	nam <- make.unique(as.character(match.call()[-1]))
	if (!missing(numeric)) nam <- nam[-length(nam)]
	if (!missing(prob)) nam <- nam[-length(nam)]
	out <- as.data.frame(table(...), stringsAsFactors=FALSE)
	for (i in 1:(ncol(out)-1))
		if (numeric) out[, i] <- as.numeric(out[, i])
	names(out)[1:(ncol(out)-1)] <- as.character(nam)
	if (prob) out$Freq <- out$Freq/sum(out$Freq)
	return(out)
}


interp <- function(x, y, z, df="BIC", raw=FALSE) {
	pred <- approx(x=x, y=y, xout=z)$y
	if (!raw) {
		offlimit <- is.na(pred)
		if (any(offlimit)) {
			pred[offlimit] <- nslm(x, y, z[offlimit], df=df)$predict
		}
	}
	return(pred)
}


cdfplot <- function(x, ...)	UseMethod("cdfplot")


cdfplot.default <- function(x, add=FALSE, ylab="F(x)", xlab="",
														lwd=1, main="", ...) {
	x <- sort(x)
	n <- length(x)
	nx <- (1:n)/n
	if (!add) {
		plot(x, nx, ylab=ylab, xlab=xlab,
				 type="s", lwd=lwd, main=main, ...)
	} else lines(x, nx, lwd=lwd, ...)
	invisible(tapply(nx, x, sum))
}

