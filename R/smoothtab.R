

smoothtab <- function(x, y, presmoothing=FALSE, postsmoothing=FALSE,
                      bandwidth=.33, grid=200, lldeg=4, llxdeg=1,
                      raw=TRUE, cdf=TRUE, margin=0.5) {
  
  if (postsmoothing && bandwidth == "auto")
    bandwidth <- dpik(x, gridsize=grid)
  
  if (missing(y) && is.vector(x)) {
    
    x <- x[!is.na(x)]
    
    ft <- conttab(x)
    colnames(ft) <- c("x", "y")
    
    if (presmoothing)
      ft[, 2] <- glm(ft[, 2] ~ poly(ft[, 1], degree=lldeg, raw=raw),
                     family=poisson)$fitted
    
    if (postsmoothing) {
      x <- rep(ft[, 1], times=round(ft[, 2]))
      
      ft <- as.data.frame(sapply(bkde(x, bandwidth=bandwidth,
                                      gridsize=grid,
                                      range.x=c(min(x)-margin, max(x)+margin),
                                      truncate=TRUE), cbind))
    }
    
    if (cdf) ft$y <- cumsum(ft[, 2])/sum(ft[, 2])
    
    out <- as.data.frame(ft)
    out <- as.smoothtab(out, presmoothing=presmoothing,
                        postsmoothing=postsmoothing,
                        bandwidth=bandwidth, grid=grid,
                        lldeg=lldeg, llxdeg=llxdeg,
                        raw=raw, cdf=cdf, margin=margin,
                        design="EG", range=c(min(x), max(x)))
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
      colnames(x.ft) <- c("x", "y")
      y.ft <- as.data.frame(cbind(sort(unique(y)),
                                  tapply(ft[, 3], ft[, 2], sum)))
      colnames(y.ft) <- c("x", "y") 
    } else {
      x.ft <- conttab(x)
      colnames(x.ft) <- c("x", "y")
      y.ft <- conttab(y)
      colnames(y.ft) <- c("x", "y")
    }
    
    if (postsmoothing) {
      if (presmoothing) {
        n <- length(x)+nrow(ft)
        prob <- round((ft[, 3]/sum(ft[, 3]))*n)
        index <- rep(1:nrow(ft), times=prob)
        xy <- ft[index, 1:2]
        colnames(xy) <- c("x", "y")
      } else xy <- cbind(x, y)
      
      ft2d <- bkde2D(xy, bandwidth=bandwidth,
                     gridsize=c(grid, grid), truncate=TRUE,
                     range.x=list(c(min(x)-margin, max(x)+margin),
                                  c(min(y)-margin, max(y)+margin)))
      
      x.ft <- data.frame(x=ft2d$x1, y=rowSums(ft2d$fhat))      
      y.ft <- data.frame(x=ft2d$x2, y=colSums(ft2d$fhat))
    }
    
    if (cdf) {
      x.ft[, 2] <- cumsum(x.ft[, 2])/sum(x.ft[, 2])
      y.ft[, 2] <- cumsum(y.ft[, 2])/sum(y.ft[, 2])
    }
    
    out <- as.smoothtab(x.ft, y.ft, presmoothing=presmoothing,
                        postsmoothing=postsmoothing,
                        bandwidth=bandwidth, grid=grid,
                        lldeg=lldeg, llxdeg=llxdeg,
                        raw=raw, cdf=cdf, margin=margin,
                        design="SG", range=c(min(x), max(x)),
                        range.y=c(min(y), max(y)))
    
    return(out)
  }
}


as.smoothtab <- function(x, y, presmoothing=FALSE, postsmoothing=FALSE,
                         bandwidth=NA, grid=NA, lldeg=NA, llxdeg=NA,
                         raw=TRUE, cdf=TRUE, margin=0.5, design,
                         range, range.y) {
  
  if (!is.data.frame(x)) {
    stop("'x' must be a data.frame.")
  } else if (ncol(x) != 2) {
    stop("Number of columns in 'x' is not 2.")
  } else names(x) <- c("x", "y")
  
  out <- list()
  
  if (missing(y)) {
    
    out$data <- x
    
    if (missing(range)) {
      out$range <- c(x[1, 1], x[nrow(x), 1])
    } else out$range <- range
    
  } else {
    
    if (!is.data.frame(y)) {
      stop("'y' must be a data.frame.")
    } else if (ncol(y) != 2) {
      stop("Number of columns in 'y' is not 2.")
    } else names(y) <- c("x", "y")
    
    out$xdata <- x
    out$ydata <- y
    
    if (missing(range)) {
      out$xlim <- c(x[1, 1], x[nrow(x), 1])
    } else out$xlim <- range
    
    if (missing(range.y)) {
      out$ylim <- c(x[1, 1], x[nrow(x), 1])
    } else out$ylim <- range.y
    
  }
  
  if (missing(design)) {
    out$design <- "none"
  } else out$design <- design
  
  class(out) <- c("smoothtab")
  out$cdf <- cdf
    
  out$presmoothing <- presmoothing
  
  if (presmoothing) {
    out$lldeg <- lldeg
    if (design == "SG") out$llxdeg <- llxdeg
    out$raw <- raw
  }
  
  out$postsmoothing <- postsmoothing
  
  if (postsmoothing) {
    out$bandwidth <- bandwidth
    out$grid <- grid
  }
  
  return(out)
}


conttab <- function(..., numeric=TRUE) {
  nam <- make.unique(as.character(match.call()[-1]))
  out <- as.data.frame(table(...), stringsAsFactors=FALSE)
  for (i in 1:(ncol(out)-1))
    if (numeric) out[, i] <- as.numeric(out[, i])
  names(out)[1:(ncol(out)-1)] <- as.character(nam)
  return(out)
}

