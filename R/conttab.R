

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

conttab <- function(..., prob=FALSE) {
  nam <- make.unique(as.character(match.call()[-1]))
  if (!missing(prob)) nam <- nam[-length(nam)]
  out <- as.data.frame(table(...), stringsAsFactors=FALSE)
  out <- as.data.frame(sapply(out, as.numeric))
  colnames(out)[1:(ncol(out)-1)] <- nam
  if (prob) out$Freq <- out$Freq/sum(out$Freq)
  return(out)
}

