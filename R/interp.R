

#' Linear interpolartion function 
#'  
#' Return a list of points which linearly interpolate given data points.
#' The missing values are estimated using the nslm function.
#' 
#' @param x,y numeric vectors giving the coordinates of the points to be interpolated.
#' @param z   an optional set of numeric values specifying where interpolation is to take place.
#' @param df  If \code{raw=FALSE} it defines the way how the number of degrees of freedom for the
#'            Natural Cubic Spline used in nslm for estimating the missing values using AIC, BIC
#'            or the sum of squared errors. Possible values are: "BIC" (default), "AIC", "LSE"
#'            (for sum of squared errors) or numeric value for df given by the user. 
#' @param raw If \code{FALSE} (default) it uses nslm for estimating the missing values, otherwise
#'            it uses approx only.
#'            
#' @seealso \code{\link{nslm}}, \code{\link{approx}}, \code{\link{ns}}
#' 
#' @examples
#' 
#' data(Tests)
#' 
#' x <- Tests[Tests$Sample == "P", "x"]
#' y <- Tests[Tests$Sample == "P", "y"]
#' 
#' st.x <- smoothtab(x)$table
#' st.y <- smoothtab(y)$table
#' 
#' interp(st.x$prob, st.x$score, st.y$prob, raw=TRUE)
#' interp(st.x$prob, st.x$score, st.y$prob)
#' 
#' @export

interp <- function(x, y, z, df="BIC", raw=FALSE) {
	pred <- approx(x=x, y=y, xout=z)$y
	if (!raw) {
		offlimit <- is.na(pred)
		if (any(offlimit)) {
			pred[offlimit] <- nslm(x, y, z[offlimit], df=df)$predict
		}
	}
	return(pred)
}

