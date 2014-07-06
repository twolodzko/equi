

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


nslm <- function(x, y, z=x, df="BIC") {
  out <- list()
  if (df %in% c("AIC", "BIC", "LSE")) {
    k <- length(unique(y))
    crit <- NULL
    for (i in 1:(k-1)) {
      fit <- lm(y~ns(x, df=i))
      if (df == "LSE") {
        pred <- predict(fit, data.frame(x=x))
        crit[i] <- sum((y-pred)^2)
      } else if (df == "AIC") { crit[i] <- AIC(fit)
      } else crit[i] <- BIC(fit)
    }
    out$choice.method <- df
    out$criteria <- crit
    out$df <- which.min(crit)
  } else {
    out$choice.method <- "none"
    out$df <- df
  }
  out$model <- lm(y~ns(x, df=out$df))
  out$predicted <- predict(out$model, data.frame(x=z))
  class(out) <- "nslm"
  return(out)
}

