

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


summary.lin <- function(x, digits=3) {
  cat(paste("y = ", round(x$intercept, digits=digits), " + ",
            round(x$slope, digits=digits), " * x", sep=""))
}


print.smoothtab <- function(x, ...) {
  if (x$design == "SG") {
    cat("1) Table for 'x'\n")
    print(x$xdata$table, ...)
    if (!is.null(x$xdata$postsmoothing))
      cat(paste("\nKernel Smoothing bandwidth:",
                round(x$xdata$postsmoothing$h, 3), "\n"))
    cat("\n2) Table for 'y'\n")
    print(x$ydata$table, ...)
    if (!is.null(x$ydata$postsmoothing))
      cat(paste("\nKernel Smoothing bandwidth:",
                round(x$ydata$postsmoothing$h, 3), "\n"))
  } else {
    print(x$table, ...)
    if (!is.null(x$postsmoothing))
      cat(paste("\nKernel Smoothing bandwidth:",
                round(x$postsmoothing$h, 3)))
  }
}


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


cdfplot.smoothtab <- function(x, add=FALSE, ...)
  plot.smoothtab(x, add=add, ...)


print.kernsmooth <- function(x, ...) {
  print(x$data, ...)
  cat(paste("\nKernel Smoothing bandwidth:",
            round(x$h, 3)))
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


print.nslm <- function(x, ...) print(x$model, ...)


summary.nslm <- function(x, ...) {
  cat(paste("Natural Cubic Spline linear regression with",
            x$df, "degrees of freedom"))
  if (x$choice.method != "none")
    cat(paste(" choosen using", x$choice.method))
  cat(".\n")
  print(x$model, ...)
  print(summary(x$model), ...)
}

