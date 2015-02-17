

#' Bootstrap standard error of equating (SEE)
#' 
#' Compute bootstrap Standard Error of Equating (SEE). 
#' 
#' @param data     a data frame.
#' @param fun      a function to pass to bootstrap. The function must return
#'                 a concordance table with source test scores in first column
#'                 named "x", and equated scores in the second column named "yx"
#'                 as in \code{\link{equi}} or \code{\link{lin}} functions.
#' @param reps     number of bootstrap repetitions (100 by default).
#' @param clusters vector or list of vectors of cluster groupping variables.
#' 
#' @return
#' 
#' Returns Bias, SE and RMSE values for each score point and averaged SEE
#' values coputed with mean or average weighted on score points probabilities. 
#'  
#'  
#' @references
#' 
#' Efron, B. & Tibshirani, R.J. (1993). An Introduction to the Bootstrap.
#' London: Chapman & Hall/CRC.
#' 
#' Field, C.A. & Welsh, A.H. (2007). Bootstrapping clustered data.
#' Journal of the Royal Statistical Society:
#' Series B (Statistical Methodology), 69(3), 369-390.
#' 
#' Kolen, M.J. & Brennan, R.J. (2004). Test Equating, Scaling, and Linking:
#' Methods and Practices. New York: Springer-Verlag.
#' 
#' Rena, S., Lai, H., Tong, W., Aminzadeh, M., Hou, X. & Lai, S. (2010).
#' Nonparametric bootstrapping for hierarchical data. Journal of
#' Applied Statistics, 37(9), 1487-1498.
#' 
#' von Davier, A.A., Holland, P.W. & Thayer, D.T. (2004). The Kernel Method
#' of Test Equating. New York: Springer-Verlag.
#' 
#' Wang, C. (2011). An investigation of bootstrap methods for estimating 
#' the standard error of equating under the common-item nonequivalent groups
#' design. doctoral PhD diss., University of Iowa. http://ir.uiowa.edu/etd/1188 
#' 
#' @seealso \code{\link{equi}}, \code{\link{smoothtab}}
#' 
#' @examples
#' 
#' data(Tests)
#' x <- Tests[Tests$Sample == "P", "x"]
#' y <- Tests[Tests$Sample == "P", "y"]
#' 
#' data <- data.frame(x=x, y=y)
#' 
#' # a function to be passed to bootstrap
#' 
#' myfun1 <- function(data) {
#'   eq <- equi(smoothtab(data$x, data$y, presmoothing=TRUE))
#'   return(eq$conc)
#' }
#' 
#' # note: the number of iterations is small
#' #       only to run faster as an example
#' 
#' see(data, myfun1, reps=25)
#' 
#' ## Bootstrap for NEAT-CE design
#' 
#' data(Tests)
#' 
#' myfun2 <- function(Tests) {
#'   p <- Tests[Tests$Sample == "P", 1:2]
#'   q <- Tests[Tests$Sample == "Q", 2:3]
#'  eq <- equi(smoothtab(p), smoothtab(q))
#'   return(eq$conc)
#' }
#' 
#' see(Tests, myfun2, reps=25)
#' 
#' @export

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
		} else {
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
	
	out <- list(bootstrapSample=cmp, estimate=est,
              see=data.frame(bias=bias, se=se, rmse=rmse))
	class(out) <- "see"
	return(out)
}


#' @export

print.see <- function(x, ...) {
	print(colMeans(x$see))
}


#' @rdname see
#' 
#' @param x                      object to plot or print.
#' @param type,col,lty,xlab,ylab arguments of the plot function.
#' @param alpha                  plot lines opacity (range [0, 1]).
#' @param \dots                  potentially further arguments
#'                               passed from other methods.
#' 
#' @export

plot.see <- function(x, type="l", col=1, lty=1, ylab="yx", xlab="x", alpha=0.5, ...) {
	col <- col2rgb(col)/255
	col <- rgb(col[1], col[2], col[3], alpha=alpha)
	matplot(t(x$bootstrapSample), type=type, col=col, lty=lty, ylab=ylab, xlab=xlab, ...)
	lines(x$estimate$yx, col="red", lwd=3)
}

