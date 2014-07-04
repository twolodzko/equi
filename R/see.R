

equi.see <- function(x, y, design="sg", reps=100, presmoothing=FALSE,
                     postsmoothing=FALSE, df="BIC", bandwidth=.33,
                     truncate=TRUE, range, range.y,
                     lldeg=4, llxdeg=1, raw=TRUE, clusters) {
  
  design <- match.arg(design, c("eg", "sg", "neat"))
  
  param.smoothtab <- list(presmoothing=presmoothing,
                          bandwidth=bandwidth, lldeg=lldeg,
                          llxdeg=llxdeg)
  param.equi <- list(df=df, truncate=truncate)
  
  if (missing(range)) {
    if (design == "neat") {
      param.equi$xlim <- c(min(x[, 2]), max(x[, 2]))
      param.equi$ylim <- c(min(y[, 2]), max(y[, 2]))
    } else param.equi$range <- c(min(y), max(y))
  }
  
  if (design == "sg") {
    data <- cbind(x, y)
    n <- nrow(data)
    k <- length(unique(x))
    
    param.smoothtab$x <- x
    param.smoothtab$y <- y
    param.equi$x <- do.call(smoothtab, param.smoothtab)
    eq <- do.call(equi, param.equi)$concordance
    yx <- eq$yx
    
    cmp <- matrix(NA, reps, k)
    
    for (i in 1:reps) {
      if (missing(clusters)) {
        index <- sample(1:n, n, replace=TRUE)
      } else {
        unq <- unique(clusters)
        nu <- length(unq)
        clsamp <- sample(unq, nu, replace=TRUE)
        index <- NULL
        for (j in clsamp) index <- c(index, which(clusters == j))
      }
      samp <- data[index, ]
      param.smoothtab$x <- samp[, 1]
      param.smoothtab$y <- samp[, 2]
      param.equi$x <- do.call(smoothtab, param.smoothtab)
      tmp <- do.call(equi, param.equi)$concordance
      rownames(tmp) <- tmp$x
      cmp[i, ] <- tmp[as.character(eq$x), "yx"]
    }
  } else if (design == "neat") {
    nx <- nrow(x)
    ny <- nrow(y)
    k <- length(unique(x[, 1]))
    
    param.smoothtab.x <- param.smoothtab
    param.smoothtab.y <- param.smoothtab
    
    param.smoothtab.x$x <- x[, 1]
    param.smoothtab.x$y <- x[, 2]
    param.equi$x <- do.call(smoothtab, param.smoothtab.x)
    param.smoothtab.y$x <- y[, 1]
    param.smoothtab.y$y <- y[, 2]
    param.equi$y <- do.call(smoothtab, param.smoothtab.y)
    eq <- do.call(equi.ce, param.equi)$concordance
    yx <- eq$yx
    
    cmp <- matrix(NA, reps, k)
    
    for (i in 1:reps) {
      if (missing(clusters)) {
        index <- sample(1:nx, nx, replace=TRUE)
      } else {
        unq <- unique(clusters[[1]])
        nu <- length(unq)
        clsamp <- sample(unq, nu, replace=TRUE)
        index <- NULL
        for (j in clsamp) index <- c(index, which(clusters[[1]] == j))
      }
      x.samp <- x[index, ]
      if (missing(clusters)) {
        index <- sample(1:ny, ny, replace=TRUE)
      } else {
        unq <- unique(clusters[[2]])
        nu <- length(unq)
        clsamp <- sample(unq, nu, replace=TRUE)
        index <- NULL
        for (j in clsamp) index <- c(index, which(clusters[[2]] == j))
      }
      y.samp <- y[index, ]
      param.smoothtab.x$x <- x.samp[, 1]
      param.smoothtab.x$y <- x.samp[, 2]
      param.equi$x <- do.call(smoothtab, param.smoothtab.x)
      param.smoothtab.y$x <- y.samp[, 1]
      param.smoothtab.y$y <- y.samp[, 2]
      param.equi$y <- do.call(smoothtab, param.smoothtab.y)
      tmp <- do.call(equi.ce, param.equi)$concordance
      rownames(tmp) <- tmp$x
      cmp[i, ] <- tmp[as.character(eq$x), "yx"]
    }
  } else {
    nx <- length(x)
    ny <- length(y)
    k <- length(unique(x))
    
    param.smoothtab$x <- x
    param.equi$x <- do.call(smoothtab, param.smoothtab)
    param.smoothtab$x <- y
    param.equi$y <- do.call(smoothtab, param.smoothtab)
    eq <- do.call(equi, param.equi)$concordance
    yx <- eq$yx
    
    cmp <- matrix(NA, reps, k)
    
    for (i in 1:reps) {
      if (missing(clusters)) {
        index <- sample(1:nx, nx, replace=TRUE)
      } else {
        unq <- unique(clusters[[1]])
        nu <- length(unq)
        clsamp <- sample(unq, nu, replace=TRUE)
        index <- NULL
        for (j in clsamp) index <- c(index, which(clusters[[1]] == j))
      }
      x.samp <- x[index]
      if (missing(clusters)) {
        index <- sample(1:ny, ny, replace=TRUE)
      } else {
        unq <- unique(clusters[[2]])
        nu <- length(unq)
        clsamp <- sample(unq, nu, replace=TRUE)
        index <- NULL
        for (j in clsamp) index <- c(index, which(clusters[[2]] == j))
      }
      y.samp <- y[index]
      param.smoothtab$x <- x.samp
      param.equi$x <- do.call(smoothtab, param.smoothtab)
      param.smoothtab$x <- y.samp
      param.equi$y <- do.call(smoothtab, param.smoothtab)
      tmp <- do.call(equi, param.equi)$concordance
      rownames(tmp) <- tmp$x
      cmp[i, ] <- tmp[as.character(eq$x), "yx"]
    }
  }
  
  bias <- colMeans(cmp, na.rm=TRUE) - yx
  se   <- apply(cmp, 2, sd, na.rm=TRUE)
  rmse <- sqrt(bias^2 + se^2)
  
  if (design == "neat") {
    w <- as.vector(table(x[, 1]))
  } else w <- as.vector(table(x))
  w <- w/sum(w)
  
  out <- list(bias=bias, se=se, rmse=rmse,
              averaged.see=data.frame(bias=mean(bias),
                                      se=mean(se),
                                      rmse=mean(rmse),
                                      wbias=sum(bias*w),
                                      wse=sum(se*w),
                                      wrmse=sum(rmse*w),
                                      row.names=""),
              bootstrap.repetitions=reps)
  print(out$averaged.see)
  class(out) <- "see"
  invisible(out)
}

