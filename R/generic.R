

print.equi <- function(x) {
  print(x$conc)
  if (x$design == "EG") {
    cat("\nEquivalent Groups (EG) design.")
  } else if (x$design == "NEAT-CE") {
    cat("\nNon-Equivalent groups with Anchor\nTest - Chain Equating (NEAT-CE) design.")
  } else cat("\nSingle Group (SG) design.")
}

plot.equi <- function(x, type="l", ident=TRUE, ...) {
  plot(x$conc, type=type, ...)
  if (ident) {
    lines(c(x$conc[1, 1]-1, x$conc[nrow(x$conc), 1]+1),
          c(x$conc[1, 2]-1, x$conc[nrow(x$conc), 2]+1), lty=2)
    legend("topleft", legend=c("ident", "equi"), lty=c(1, 2),
           bty="n", cex=0.8)
  }
}

summary.lin <- function(x) {
  cat(paste("y = ", round(x$intercept, 3), " + ",
            round(x$slope, 3), " * x", sep=""))
}

print.smoothtab <- function(x) {
  if (x$design == "SG") {
    cat("1) Table for 'x'\n")
    print(x$xdata$table)
    if (!is.null(x$xdata$postsmoothing))
      cat(paste("\nKernel Smoothing bandwidth:",
                round(x$xdata$postsmoothing$h, 3), "\n"))
    cat("\n2) Table for 'y'\n")
    print(x$ydata$table)
    if (!is.null(x$ydata$postsmoothing))
      cat(paste("\nKernel Smoothing bandwidth:",
                round(x$ydata$postsmoothing$h, 3), "\n"))
  } else {
    print(x$table)
    if (!is.null(x$postsmoothing))
      cat(paste("\nKernel Smoothing bandwidth:",
                round(x$postsmoothing$h, 3)))
  }
}

plot.smoothtab <- function(x, type="l", ...) {
  if (x$design == "SG") {
    lim <- c(min(x$xdata$table$score, x$ydata$table$score),
             max(x$xdata$table$score, x$ydata$table$score))
    plot(x$xdata$table, type=type,
         xlab="", xlim=lim, ...)
    lines(x$ydata$table, type=type, lty=2, ...)
  } else plot(x$table, type=type, xlab="", ...)
}

cdfplot.smoothtab <- function(x, ...) plot(x, ...)

print.kernsmooth <- function(x) {
  print(x$data)
  cat(paste("\nKernel Smoothing bandwidth:",
            round(x$h, 3)))
}

plot.kernsmooth <- function(x, type="l", ...) {
  plot(x$data, type=type, ...)
}

cdfplot.kernsmooth <- function(x, type="l", ...) {
  plot(x$data$x, cumsum(x$data$Freq)/sum(x$data$Freq), 
       type=type, ylab="F(x)", xlab="")
}

print.nslm <- function(x) print(x$model)

summary.nslm <- function(x) {
  cat(paste("Natural Cubic Spline linear regression with",
            x$df, "degrees of freedom"))
  if (x$choice.method != "none")
    cat(paste(" choosen using", x$choice.method))
  cat(".\n")
  print(x$model)
  print(summary(x$model))
}

