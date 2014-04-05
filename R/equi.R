

equi <- function(x, y, range, fun="none", df=10,
                 truncate=TRUE, margin=0.5) {

  fun <- match.arg(fun, c("none", "ns", "bs"))

  if (!missing(y) && class(y) == "equi") {
    return(equi.predict(x, y, truncate=truncate,
                        fun=y$fun, range=range, df=df))  
  } else {
    
    if (missing(y) && class(x) == "smoothtab"
        && x$design == "SG") {
      
      xdata <- x$xdata
      ydata <- x$ydata
      if (missing(range)) range <- x$ylim
      design <- "SG"
      
    } else if (class(x) == "smoothtab"
               && class(y) == "smoothtab"
               && x$design == "SG"
               && y$design == "SG") {
      
      if (missing(range)) range <- y$ylim
      
      return(equi.ce(x, y, ylim=range, fun=fun,
                     df=df, truncate=truncate))
      
    } else if (class(x) == "smoothtab"
               && class(y) == "smoothtab"
               && x$design == "EG"
               && y$design == "EG") {
      
      xdata <- x$data
      ydata <- y$data
      if (missing(range)) range <- y$range
      design <- "EG"
      
    } else stop("'x' and 'y' are supposed to be 'smoothtab' objects.")
    
    if (fun %in% c("ns", "bs")) {
      if (fun == "ns") {
        fx <- lm(y~ns(x, df=min(df, length(unique(x))-1)), data=xdata)
        fy <- lm(x~ns(y, df=min(df, length(unique(y))-1)), data=ydata)
      } else {
        fx <- lm(y~bs(x, df=min(df, length(unique(x))-1)), data=xdata)
        fy <- lm(x~bs(y, df=min(df, length(unique(y))-1)), data=ydata)
      }
      x.ft <- predict(fx, data.frame(x=xdata$x))
      pred <- predict(fy, data.frame(y=x.ft))
      conc <- data.frame(x=xdata$x, yx=pred)
    } else {
      pred <- approx(x=ydata$y, y=ydata$x, xdata$y)$y
      conc <- data.frame(x=xdata$x, yx=pred)
    }
    
    if (truncate)
      conc$yx <- trun(conc$yx, range=range, margin=margin)
    
    out <- list(concordance=conc, fun=fun,
                range=range, truncate=truncate,
                design=design)
    
    out$presmoothing <- x$presmoothing
    out$postsmoothing <- x$postsmoothing
    
    if (x$presmoothing) {
      out$lldeg <- x$lldeg
      out$llxdeg <- x$llxdeg
      out$raw <- x$raw
    }
    
    if (x$postsmoothing) {
      out$bandwidth <- x$bandwidth
      out$grid <- x$grid
    }
    
    if (fun %in% c("ns", "bs")) out$df <- df
    class(out) <- "equi"
    return(out)
  }
}


equi.ce <- function(x, y, truncate=TRUE, fun="none",
                    df=10, xlim, ylim, margin=0.5) {
  
  fun <- match.arg(fun, c("none", "ns", "bs"))

  if (missing(xlim)) xlim <- x$ylim
  if (missing(ylim)) ylim <- y$ylim
  
  range <- ylim
  
  fx <- equi(x, fun=fun, truncate=truncate,
             range=xlim, df=df, margin=margin) 
  fy <- equi(y, fun=fun, truncate=truncate,
             range=ylim, df=df, margin=margin)
  yx <- equi(equi(x$xdata$x, fx), fy)
  
  out <- list(concordance=data.frame(x=x$xdata$x, yx=yx),
              fun=fun, range=range, truncate=truncate,
              design="NEAT-CE")
  
  out$presmoothing <- x$presmoothing
  out$postsmoothing <- x$postsmoothing
  
  if (x$presmoothing) {
    out$lldeg <- x$lldeg
    out$llxdeg <- x$llxdeg
    out$raw <- x$raw
  }
  
  if (x$postsmoothing) {
    out$bandwidth <- x$bandwidth
    out$grid <- x$grid
  }
  
  if (fun %in% c("ns", "bs")) out$df <- df
  class(out) <- "equi"
  return(out)
}


equi.predict <- function(x, y, truncate, fun, range, df=10, margin=0.5) {
  
  if (missing(fun)) fun <- y$fun
  
  if (fun == "ns") {
    yx <- as.numeric(predict(lm(yx~ns(x, df=min(df, length(unique(x))-1)),
                                data=y$concordance), data.frame(x=x)))
  } else if (fun == "bs") {
    yx <- as.numeric(predict(lm(yx~ns(x, df=min(df, length(unique(x))-1)),
                                data=y$concordance), data.frame(x=x)))
  } else yx <- approx(x=y$concordance$x, y=y$concordance$yx, x)$y
  
  if (missing(truncate)) truncate <- y$truncate
  if (truncate) yx <- trun(yx, y$range, margin=margin)
  
  return(yx)
}

