

smoothtab <- function(x, y, presmoothing=FALSE, postsmoothing=FALSE,
                      bandwidth="auto", lldeg=4, llxdeg=1,
                      raw=TRUE, cdf=TRUE, margin=0.5) {
  
  if (missing(y) && is.vector(x)) {
    
    x <- x[!is.na(x)]
    
    ft <- conttab(x)
    colnames(ft) <- c("x", "y")
    
    if (presmoothing)
      ft[, 2] <- glm(ft[, 2] ~ poly(ft[, 1], degree=lldeg, raw=raw),
                     family=poisson)$fitted
    
    if (postsmoothing) {
      ft[, 2] <- ft[, 2]/sum(ft[, 2])
      ks <- kern(x, xj=ft[, 1], p=ft[, 2], m=mean(x), s=var(x), h=bandwidth)
      bandwidth <- ks$h
      ft[, 2] <- ks$data$Freq
    }
    
    if (cdf) ft$y <- cumsum(ft[, 2])/sum(ft[, 2])
    
    out <- as.data.frame(ft)
    out <- as.smoothtab(out, presmoothing=presmoothing,
                        postsmoothing=postsmoothing,
                        bandwidth=bandwidth, 
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
      x.ft[, 2] <- x.ft[, 2]/sum(x.ft[, 2])
      ks <- kern(x, xj=x.ft[, 1], p=x.ft[, 2], m=mean(x), s=var(x), h=bandwidth)
      x.bandwidth <- ks$h
      x.ft[, 2] <- ks$data$Freq
      
      y.ft[, 2] <- y.ft[, 2]/sum(y.ft[, 2])
      ks <- kern(y, xj=y.ft[, 1], p=y.ft[, 2], m=mean(y), s=var(y), h=bandwidth)
      y.bandwidth <- ks$h
      y.ft[, 2] <- ks$data$Freq
    }
    
    if (cdf) {
      x.ft[, 2] <- cumsum(x.ft[, 2])/sum(x.ft[, 2])
      y.ft[, 2] <- cumsum(y.ft[, 2])/sum(y.ft[, 2])
    }
    
    out <- as.smoothtab(x.ft, y.ft, presmoothing=presmoothing,
                        postsmoothing=postsmoothing,
                        bandwidth=x.bandwidth,
                        y.bandwidth=y.bandwidth, 
                        lldeg=lldeg, llxdeg=llxdeg,
                        raw=raw, cdf=cdf, margin=margin,
                        design="SG", range=c(min(x), max(x)),
                        range.y=c(min(y), max(y)))
    
    return(out)
  }
}


as.smoothtab <- function(x, y, presmoothing=FALSE, postsmoothing=FALSE,
                         bandwidth=NA, y.bandwidth, lldeg=NA, llxdeg=NA,
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
    if (!missing(y)) {
      out$x.bandwidth <- bandwidth
      out$y.bandwidth <- ifelse(missing(y.bandwidth), bandwidth, y.bandwidth)
    } else out$bandwidth <- bandwidth
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


kern <- function(x, xj, p, m, s, h="auto", hmin=0.1, hmax=1, pen=FALSE, grid=0.1) {
  if (h == "auto") {
    f <- function(h) kern(x, xj, p, m, s, h, pen=TRUE)$pen
    return(kern(x, xj, p, m, s, h=optimize(f, lower=hmin, upper=hmax)$minimum))
  } else {
    k <- length(xj)
    a <- sqrt(s/(s+h^2))
    
    res <- numeric(k)
    for (i in 1:k) {
      Rjx <- (x - a*xj[i] - (1-a)*m) / (a*h)
      res[i] <- sum(p * dnorm(Rjx) / (a*h))
    }
    
    out <- list(data=data.frame(x=xj, Freq=res), h=h)
    class(out) <- "kernsmooth"
    
    if (pen) {
      der <- numeric(k)
      for (i in 1:k) {
        Rjx <- (x - a*xj[i] - (1-a)*m) / (a*h)
        der[i] <- -sum(p * dnorm(Rjx) / (a*h)^2 * Rjx )
      }
      p1 <- sum((p-res)^2)
      p2 <- sum( (der > 0) * (1 - (der <= 0) ))
      out$pen <- p1+p2
    }
    return(out)
  }
}

