

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
      conc$yx <- trun(conc$yx, range=range, margin=margin)
    
    out <- list(concordance=conc, range=range,
                truncate=truncate, df=df, design=design)
    
    class(out) <- "equi"
    return(out)
  }
  UseMethod("equi")
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
  if (truncate) yx <- trun(yx, y$range, margin=margin)
  return(yx)
}


print.equi <- function(x, ...) {
	print(x$conc, ...)
	if (x$design == "EG") {
		cat("\nEquivalent Groups (EG) design.")
	} else if (x$design == "NEAT-CE") {
		cat("\nNon-Equivalent groups with Anchor\nTest - Chain Equating (NEAT-CE) design.")
	} else cat("\nSingle Group (SG) design.")
}


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

