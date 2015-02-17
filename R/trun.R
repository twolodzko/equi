

#' Truncate
#' 
#' Truncate a vector given lower and upper margins. 
#' 
#' @param x           numeric vector.
#' @param lower,upper lower and upper margins.
#' @param margin      value of margin is subtracted form
#'                    lower value of range and added to
#'                    upper value of range.
#'               
#' @examples
#' 
#' trun(1:10, 4, 8, margin=0)
#' trun(1:10, 4, 8, margin=0.5)
#' 
#' @export

trun <- function(x, lower, upper, margin=0) {
  if (!missing(lower)) {
    lower <- lower - margin
    x[x < lower] <- lower
  }
  if (!missing(upper)) {
    upper <- upper + margin
    x[x > upper] <- upper
  }
  return(x)
}

