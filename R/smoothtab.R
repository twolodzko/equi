

#' Presmoothing and postsmoothing of empirical distributions
#' 
#' Function for Log-linear presmoothing and/or Gaussian Kernel postsmoothing. 
#' 
#' @param x, y          numeric vectors.
#' @param presmoothing  if \code{TRUE} Log-linear presmoothing is applied.
#' @param postsmoothing if \code{TRUE} Gaussian Kernel postsmoothing is applied.
#' @param bandwidth     sets bandwidth for Kernel Smoothing. Use "auto"
#'                      (default) for automatic selection of bandwidth. 
#' @param lldeg         degree of the polynomial in log-linear presmoothing
#'                      (deafult is 4).
#' @param llxdeg        degree of the polynomial in log-linear presmoothing
#'                      for interaction term (deafult is 1).
#' @param raw           if \code{TRUE} computes raw polynomials for log-linear
#'                      presmoothing, see: \code{\link{poly}}.
#' @param cdf           if \code{FALSE} compute probabilities rather than cumulative
#'                      probabilities (default).
#' @param margin        if \code{postsmoothing=TRUE}, it defines the margins of points
#'                      range to be created
#'                      using Gaussian Kernel postsmoothing. The function returns
#'                      point values +/- margin.
#' @param grid          if \code{postsmoothing=TRUE}, it defines the number of points to
#'                      be created using Gaussian Kernel postsmoothing.
#'
#' @return
#' 
#' Returns two-column \code{data.frame} with unique score points and coresponding probabilities,
#' or a list of two \code{data.frames} for marginal probabilities of joint distributions. 
#' 
#' @note
#' See also \pkg{equate} and \pkg{kequate} packages.
#' 
#' @references
#' 
#' Holland, P.W. & Thayer, D.T. (2000). \emph{Univariate and Bivariate Loglinear Models
#' for Discrete Test Score Distributions.} Journal of Educational and Behavioral
#' Statistics, 25(2), 133-183.
#' 
#' Kolen, M.J. & Brennan, R.J. (2004). \emph{Test Equating, Scaling, and Linking:
#' Methods and Practices.} New York: Springer-Verlag.
#' 
#' von Davier, A.A., Holland, P.W. & Thayer, D.T. (2004). \emph{The Kernel Method
#' of Test Equating.} New York: Springer-Verlag.
#' 
#' Wand, M.P. & Jones, M.C. (1995). \emph{Kernel Smoothing.} London: Chapman & Hall/CRC. 
#'
#' @examples
#' 
#' data(Tests)
#' x <- Tests[Tests$Sample == "P", "x"]
#' y <- Tests[Tests$Sample == "P", "y"]
#' 
#' (st <- smoothtab(x))
#' plot(st)
#' 
#' smoothtab(x, presmoothing=TRUE)
#' smoothtab(x, y, presmoothing=TRUE, postsmoothing=TRUE)
#' 
#' @export

smoothtab <- function(x, y, presmoothing=FALSE, postsmoothing=FALSE,
                      bandwidth="auto", lldeg=4, llxdeg=1,
                      raw=TRUE, cdf=TRUE, margin=0.5, grid=100) {
  
  if (missing(y) && is.vector(x)) {
    
    x <- x[!is.na(x)]
    ft <- conttab(x)
    
    if (presmoothing)
      ft[, 2] <- glm(ft[, 2] ~ poly(ft[, 1], degree=lldeg, raw=raw),
                     family=poisson)$fitted
    
    if (postsmoothing) {
      ft[, 2] <- ft[, 2]/sum(ft[, 2])
      x.new <- seq(min(ft[, 1]), max(ft[, 1]), length.out=grid)
      ks <- kern(x.new, xj=ft[, 1], p=ft[, 2], m=mean(x), s=var(x), h=bandwidth)
      x.bandwidth <- ks$h
      ft <- ks$data
    }
    
    if (cdf) ft[, 2] <- cumsum(ft[, 2])/sum(ft[, 2])
    
    colnames(ft) <- c("score", "prob")
    
    out <- list(design="EG", table=as.data.frame(ft),
                range=c(min(x), max(x)), cdf=cdf)
    
    if (presmoothing) {
      out$presmoothing <- list(lldeg=lldeg, llxdeg=llxdeg, raw=raw)
    }     
    if (postsmoothing) {
      out$postsmoothing <- list(bandwidth=bandwidth,
                                h=ifelse(bandwidth=="auto",
                                         x.bandwidth, bandwidth))
    } 
    
    class(out) <- "smoothtab"
    return(out)
    
  } else {
    
    if (missing(y)) {
      y <- x[, 2]
      x <- x[, 1]
    }
    
    filter <- !is.na(x) & !is.na(y)
    x <- x[filter]
    y <- y[filter]
    
    if (presmoothing) {
      ft <- conttab(x, y)
      ft[, 3] <- glm(ft[, 3] ~ poly(ft[, 1], degree=lldeg, raw=raw)
                     + poly(ft[, 2], degree=lldeg, raw=raw)
                     + poly(ft[, 1]*ft[, 2], degree=llxdeg,
                            raw=raw), family=poisson)$fitted  
      
      x.ft <- as.data.frame(cbind(sort(unique(x)),
                                  tapply(ft[, 3], ft[, 1], sum)))
      y.ft <- as.data.frame(cbind(sort(unique(y)),
                                  tapply(ft[, 3], ft[, 2], sum)))
    } else {
      x.ft <- conttab(x)
      y.ft <- conttab(y)
    }
    
    if (postsmoothing) {
      x.ft[, 2] <- x.ft[, 2]/sum(x.ft[, 2])
      x.new <- seq(min(x.ft[, 1])-margin, max(x.ft[, 1])+margin, length.out=grid)
      ks <- kern(x.new, xj=x.ft[, 1], p=x.ft[, 2], m=mean(x), s=var(x), h=bandwidth)
      x.bandwidth <- ks$h
      x.ft <- ks$data      
      y.ft[, 2] <- y.ft[, 2]/sum(y.ft[, 2])
      y.new <- seq(min(y.ft[, 1]), max(y.ft[, 1]), length.out=grid)
      ks <- kern(y.new, xj=y.ft[, 1], p=y.ft[, 2], m=mean(y), s=var(y), h=bandwidth)
      y.bandwidth <- ks$h
      y.ft <- ks$data
    }
    
    if (cdf) {
      x.ft[, 2] <- cumsum(x.ft[, 2])/sum(x.ft[, 2])
      y.ft[, 2] <- cumsum(y.ft[, 2])/sum(y.ft[, 2])
    }
    
    colnames(x.ft) <- c("score", "prob")
    colnames(y.ft) <- c("score", "prob")
    
    out <- list()    
    out$design <- "SG"
    out$xdata <- list(table=as.data.frame(x.ft),
                      range=c(min(x), max(x)), cdf=cdf)
    
    if (presmoothing) {
      out$xdata$presmoothing <- list(lldeg=lldeg, llxdeg=llxdeg, raw=raw)
    } 
    if (postsmoothing) {
      out$xdata$postsmoothing <- list(bandwidth=bandwidth,
                                      h=ifelse(bandwidth=="auto",
                                               x.bandwidth, bandwidth))
    } 
    
    out$ydata <- list(table=as.data.frame(y.ft),
                      range=c(min(y), max(y)), cdf=cdf)
    
    if (presmoothing) {
      out$ydata$presmoothing <- list(lldeg=lldeg, llxdeg=llxdeg, raw=raw)
    } 
    if (postsmoothing) {
      out$ydata$postsmoothing <- list(bandwidth=bandwidth,
                                      h=ifelse(bandwidth=="auto",
                                               y.bandwidth, bandwidth))
    } 
    
    class(out) <- "smoothtab"
    return(out)
  }
}


#' @rdname smoothtab
#' @export

as.smoothtab <- function(x) {
  if (is.matrix(x) || is.data.frame(x)) {
    if (ncol(x) != 2) stop("Number of columns is different than 2.")
    out <- list()
    out$table <- as.data.frame(x)
    colnames(out$table) <- c("score", "prob")
    class(out) <- "smoothtab"
    return(out)
  } else stop("'x' must be a matrix or data.frame.")
}


#' @rdname smoothtab
#' @export

is.smoothtab <- function(x) {
	inherits(x, "smoothtab")
}


#' @export

print.smoothtab <- function(x, ...) {
	if (x$design == "SG") {
		cat("1) Table for 'x'\n")
		print(x$xdata$table, ...)
		if (!is.null(x$xdata$postsmoothing))
			cat(paste("\nKernel Smoothing bandwidth:",
                round(x$xdata$postsmoothing$h, digits=getOption("digits")), "\n"))
		cat("\n2) Table for 'y'\n")
		print(x$ydata$table, ...)
		if (!is.null(x$ydata$postsmoothing))
			cat(paste("\nKernel Smoothing bandwidth:",
                round(x$ydata$postsmoothing$h, digits=getOption("digits")), "\n"))
	} else {
		print(x$table, ...)
		if (!is.null(x$postsmoothing))
			cat(paste("\nKernel Smoothing bandwidth:",
                round(x$postsmoothing$h, digits=getOption("digits"))))
	}
}


#' @rdname smoothtab
#' 
#' @param type    type of the plot.
#' @param add     add a plot to previous one.
#' @param lty     a vector of line types, see: \code{\link{par}}.
#' @param \dots   potentially further arguments passed from other methods.
#' 
#' @export

plot.smoothtab <- function(x, type="s", lty=1:6, add=FALSE, ...) {
	if (x$design == "SG") {
		if (!add) {
			lim <- c(min(x$xdata$table$score, x$ydata$table$score),
               max(x$xdata$table$score, x$ydata$table$score))
			plot(x$xdata$table, type=type, xlab="",
           xlim=lim, lty=lty[1], ...)
			lines(x$ydata$table, type=type, lty=lty[2], ...)
		} else {
			lines(x$ydata$table, type=type, lty=lty[1], ...)
			lines(x$ydata$table, type=type, lty=lty[2], ...)
		}
	} else {
		if (!add) {
			plot(x$table, type=type, xlab="", lty=lty[1], ...)
		} else lines(x$table, type=type, xlab="", lty=lty[1], ...)
	}
}


#' @rdname smoothtab
#' @export

cdfplot.smoothtab <- function(x, add=FALSE, ...) {
	plot.smoothtab(x, add=add, ...)
}

