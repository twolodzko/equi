

#' Recode variable 
#' 
#' Recode a vector from chosen values to new values. Any value not existing in from vector
#' will be recoded to other. Numeric vector can be truncated. 
#' 
#' @param x         vector.
#' @param from      source values.
#' @param to        target values.
#' @param other     if provided, other values are recoded to other.
#' @param truncate  if \code{TRUE} values below min and above max of the to vector are recoded
#'                  to min and max respectively. See: \code{\link{trun}}.
#' @param margin    see: \code{\link{trun}}.
#' 
#' @seealso \code{\link{trun}}
#' 
#' @examples
#' 
#' recode(1:10, from=3:5, to=5:3, truncate=TRUE)
#' recode(1:10, from=3:5, to=5:3, other=NA)
#' 
#' @export

recode <- function(x, from, to, other, truncate=TRUE, margin=0.5) {
  new_x <- x
  k <- length(from)
  range <- c(min(to, na.rm=TRUE), max(to, na.rm=TRUE))
  for (i in 1:k) new_x[x == from[i]] <- to[i]
  if (!missing(other)) new_x[!(x %in% from)] <- other
  if (truncate) new_x <- trun(new_x, lower=range[1], upper=range[2], margin=margin)
  return(new_x)
}

