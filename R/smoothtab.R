

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


is.smoothtab <- function(x) "smoothtab" %in% class(x)


conttab <- function(..., numeric=TRUE, prob=FALSE) {
  nam <- make.unique(as.character(match.call()[-1]))
  if (!missing(numeric)) nam <- nam[-length(nam)]
  if (!missing(prob)) nam <- nam[-length(nam)]
  out <- as.data.frame(table(...), stringsAsFactors=FALSE)
  for (i in 1:(ncol(out)-1))
    if (numeric) out[, i] <- as.numeric(out[, i])
  names(out)[1:(ncol(out)-1)] <- as.character(nam)
  if (prob) out$Freq <- out$Freq/sum(out$Freq)
  return(out)
}


kern <- function(x, xj, p, m, s, h="auto", hmin=0.1, hmax=1, k=1, pen=FALSE) {
  if (h == "auto") {
    f <- function(h) kern(x, xj, p, m, s, h, k=k, pen=TRUE)$pen
    return(kern(x, xj, p, m, s, h=optimize(f, lower=hmin, upper=hmax)$minimum))
  } else {
    k <- length(xj)
    n <- length(x)
    a <- sqrt(s/(s+h^2))
    res <- NULL
    
    Rjx <- function(x, a, xj, m, h, i)
      (x[i] - a*xj - (1-a)*m) / (a*h)
    
    for (i in 1:n) {
      res[i] <- sum(p * dnorm(Rjx(x, a, xj, m, h, i)) / (a*h))
    }
    
    out <- list(data=data.frame(x=x, Freq=res), h=h)
    class(out) <- "kernsmooth"
    
    if (pen) {
      p.hat <- NULL
      der <- NULL
      for (i in 1:k) {
        p.hat[i] <- sum(p * dnorm(Rjx(xj, a, xj, m, h, i)) / (a*h))
        der[i] <- -sum(p * dnorm(Rjx(x, a, xj, m, h, i)) / (a*h) )
      }
      p1 <- sum((p-p.hat)^2)
      p2 <- sum( (der > 0) * (1 - (der <= 0) ))
      out$pen <- p1+k*p2
    } 
    return(out)
  }
}


kernsmooth <- function(x, h="auto", grid=100) {
  x <- x[!is.na(x)]
  m <- mean(x)
  s <- sd(x)
  ct <- conttab(x)
  xj <- ct$x
  p <- ct$Freq/sum(ct$Freq)
  x.new <- seq(min(xj), max(xj), length.out=grid)
  return(kern(x.new, xj, p, m, s, h=h))
}

