

#' Equipercentile equating 
#' 
#' Functions for equipercentile equating in Equivalent Groups (EG), Single Group (SG)
#' and Nonequivalent groups with Anchor Test using Chained Equating (NEAT-CE) designs. 
#' 
#' @param x        smoothtab object or numeric vector. For NEAT first column of smoothtab
#'                 data frame is source test and second is anchor test (see Details).
#' @param y        smoothtab or equi object. For NEAT first column of smoothtab data frame
#'                 is anchor test and second is target test (see Details).
#' @param range    two-item vector of minimum and maximum score points of \code{y}.
#' @param truncate if \code{TRUE} truncates the values using trun.
#' @param df       degrees of freedom for \code{ns} splines. If \code{df} is set to "BIC", "AIC" or "LSE"
#'                 it uses BIC, AIC or sum of squared errors for choosing the best df value. 
#' @param margin   see: \code{\link{trun}}.
#' 
#' @details
#' 
#' If \code{x} is smoothtab object and \code{y} is not provided, single group equating is applied.
#' If \code{x} and \code{y} are smoothtab objects, equivalent groups or NEAT-CE equating is applied,
#' depending on the kind of data provided. If \code{x} is numeric vector and \code{y} is equi object,
#' \code{x} is transformed using equating function provided in \code{y}. 
#' 
#' @return
#' 
#' Returns object of class equi, that can be used for transforming test scores using
#' equi function (see examples). 
#' 
#' @references
#' 
#' Kolen, M.J. & Brennan, R.J. (2004). \emph{Test Equating, Scaling, and Linking:
#' Methods and Practices.} New York: Springer-Verlag.
#' 
#' von Davier, A.A., Holland, P.W. & Thayer, D.T. (2004). \emph{The Kernel Method
#' of Test Equating.} New York: Springer-Verlag.
#' 
#' Green, P.J. & Silverman, B.W. (1993). \emph{Nonparametric Regression and Generalized
#' Linear Models: A roughness penalty approach.} London: Chapman & Hall/CRC. 
#' 
#' @seealso \code{\link{ns}}, \code{\link{bs}}, \code{\link{poly}}, \code{\link{approx}}
#' 
#' @examples
#' 
#' # Single Group design with presmoothing
#' 
#' data(Tests)
#' x <- Tests[Tests$Sample == "P", "x"]
#' y <- Tests[Tests$Sample == "P", "y"]
#' 
#' (eq <- equi(smoothtab(x, y, presmoothing=TRUE)))
#' yx <- equi(x, eq)
#' 
#' moments(list(x, y, yx))
#' 
#' # Equivalent Groups design
#' 
#' data(Tests)
#' x <- Tests[Tests$Sample == "P", "x"]
#' y <- Tests[Tests$Sample == "P", "y"]
#' 
#' (eq <- equi(smoothtab(x), smoothtab(y)))
#' yx <- equi(x, eq)
#' 
#' # NEAT-CE design, 'x' is equated to 'y' via anchor test 'a' (xa, ya).
#' # in 'p' first column is source test, second is anchor test
#' # in 'q' first column is anchor test, second is target test
#' 
#' data(Tests)
#' p <- Tests[Tests$Sample == "P", 1:2]
#' q <- Tests[Tests$Sample == "Q", 2:3]
#' 
#' (eq <- equi(smoothtab(p), smoothtab(q)))
#' yx <- equi(p[, 1], eq)
#' 
#' @export

equi <- function(x, y, range, truncate=TRUE, df="BIC", margin=0.5) {

  if (!missing(y) && class(y) == "equi") {
    return(equi.predict(x, y, truncate=truncate,
                        range=range, df=df))  
  } else {
    
    if (missing(y) && class(x) == "smoothtab"
        && x$design == "SG") {
      
      xdata <- x$xdata$table
      ydata <- x$ydata$table
      if (missing(range)) range <- x$ydata$range
      design <- "SG"
      
    } else if (class(x) == "smoothtab"
               && class(y) == "smoothtab"
               && x$design == "SG"
               && y$design == "SG") {
      
      if (missing(range)) range <- y$ydata$range
      
      return(equi.ce(x, y, ylim=range, df=df, truncate=truncate))
      
    } else if (class(x) == "smoothtab"
               && class(y) == "smoothtab"
               && x$design == "EG"
               && y$design == "EG") {
      
      xdata <- x$table
      ydata <- y$table
      if (missing(range)) range <- y$range
      design <- "EG"
      
    } else stop("'x' and 'y' are supposed to be 'smoothtab' objects.")
    
    conc <- data.frame(x=xdata$score, yx=interp(ydata$prob, ydata$score, xdata$prob, df=df))
    
    if (truncate)
      conc$yx <- trun(conc$yx, lower=range[1], upper=range[2], margin=margin)
    
    out <- list(concordance=conc, range=range,
                truncate=truncate, df=df, design=design)
    
    class(out) <- "equi"
    return(out)
  }
}


equi.ce <- function(x, y, truncate=TRUE, df="BIC", xlim, ylim, margin=0.5) {
  
  if (missing(xlim)) xlim <- x$ylim
  if (missing(ylim)) ylim <- y$ylim
  
  range <- ylim
  
  fx <- equi(x, truncate=truncate,
             range=xlim, df=df, margin=margin) 
  fy <- equi(y, truncate=truncate,
             range=ylim, df=df, margin=margin)
  yx <- equi(equi(x$xdata$table$score, fx), fy)
  
  out <- list(concordance=data.frame(x=x$xdata$table$score, yx=yx),
              range=range, truncate=truncate, df=df, design="NEAT-CE")
  
  class(out) <- "equi"
  return(out)
}


equi.predict <- function(x, y, truncate, range, df="BIC", margin=0.5) {
  yx <- interp(y$concordance$x, y$concordance$yx, x, df=df)
  if (missing(truncate)) truncate <- y$truncate
  if (truncate) yx <- trun(yx, lower=y$range[1], upper=y$range[2], margin=margin)
  return(yx)
}


#' @export

print.equi <- function(x, ...) {
	print(x$conc, ...)
	if (x$design == "EG") {
		cat("\nEquivalent Groups (EG) design.")
	} else if (x$design == "NEAT-CE") {
		cat("\nNon-Equivalent groups with Anchor\nTest - Chain Equating (NEAT-CE) design.")
	} else cat("\nSingle Group (SG) design.")
}


#' @rdname equi 
#' 
#' @param diff  plot a difference from identity function.
#' @param ref   plot a reference line, i.e. an identity function.
#' @param type  type of the plot.
#' @param \dots potentially further arguments passed from other methods.
#' 
#' @export

plot.equi <- function(x, diff=FALSE, ref=TRUE, type="l", ...) {
	if (!diff) {
		plot(x$conc, type=type, ...)
		if (ref)
			lines(c(x$conc[1, 1]-1, x$conc[nrow(x$conc), 1]+1),
						c(x$conc[1, 2]-1, x$conc[nrow(x$conc), 2]+1), lty=2)
	} else {
		plot(x$conc$x, x$conc$yx-x$conc$x, type=type,
				 xlab="x", ylab="yx - x", ...)
		title("Difference from identity function")
		if (ref) abline(h=0, lty=2)
	}
}

