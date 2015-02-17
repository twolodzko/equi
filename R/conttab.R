

#' Create contingency table 
#' 
#' Create a contingency table from given variables. 
#' 
#' @param \dots    contingency table is build of variables in \dots
#' @param prob     returns probabilities instead of counts.
#' 
#' @seealso \code{\link{table}}
#' 
#' @examples
#' 
#' data(Tests)
#' x <- Tests[Tests$Sample == "P", "x"]
#' y <- Tests[Tests$Sample == "P", "y"]
#' 
#' conttab(x)
#' conttab(x, y)
#' 
#' @export

conttab <- function(..., prob=FALSE, unique=FALSE) {
  nam <- make.unique(as.character(match.call()[-1]))
  if (!missing(prob)) nam <- nam[-length(nam)]
  if (!missing(unique)) nam <- nam[-length(nam)]
  data <- as.matrix(data.frame(...))
  out <- .Call('equi_count_cpp', PACKAGE = 'equi', data, unique=FALSE, freq=!prob)
  out <- as.data.frame(out)
  colnames(out) <- c(nam, "Freq")
  return(out)
}


# conttab.pureR <- function(..., numeric=TRUE, prob=FALSE) {
#   nam <- make.unique(as.character(match.call()[-1]))
#   if (!missing(numeric)) nam <- nam[-length(nam)]
#   if (!missing(prob)) nam <- nam[-length(nam)]
#   out <- as.data.frame(table(...), stringsAsFactors=FALSE)
#   for (i in 1:(ncol(out)-1))
#       if (numeric) out[, i] <- as.numeric(out[, i])
#   names(out)[1:(ncol(out)-1)] <- as.character(nam)
#   if (prob) out$Freq <- out$Freq/sum(out$Freq)
#   return(out)
# }

