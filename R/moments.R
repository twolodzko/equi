

#' Moments
#' 
#' Calculates moments of given distribution up to the nth moment.
#' 
#' @param x      numeric vector or list of vectors.
#' @param num    calculate up to the nth moment given by \code{num}.
#' @param excess if \code{TRUE} calculates excess kurtosis.
#' 
#' @examples
#' 
#' moments(rpois(1000, 5))
#' moments(rnorm(1000), num=6)
#' 
#' @export

moments <- function(x, num=4, excess=TRUE) {
  if (is.list(x))
    return(sapply(x, moments, num=num, excess=excess))
  x <- x[!is.na(x)]
  m <- mean(x)
  s <- sd(x)
  out <- array(c(m, s))
  if (num > 2) out[3] <- mean(((x-m)/s)^3)
  if (num > 3) {
    out[4] <- mean(((x-m)/s)^4)
    if (excess) out[4] <- out[4] - 3
  }
  if (num > 4) for (i in 5:num) out[i] <- mean((x-m)^i)
  names(out) <- c("mean", "sd", "skew", "kurt", 5:num)[1:num] 
  return(out)
}

