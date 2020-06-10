#' Testing whether a given year is out-of-control
#' 
#' This function performs hypothesis testing for monitoring whether a given year is
#' out-of-control or in-control
#' 
#' @author Daeyoung Lim, \email{daeyoung.lim@uconn.edu}
#' @param y_surveill a vector of measurements to monitor
#' @param yh_bar a vector of historical measurements to estimate the parameters
#' @param sigma the estimate of the standard deviation
#' @param nyears # of years for the historical data
#' @param sig_level significance level
#' @param C the constant C that determines the power
#' @param B1 # of outer loop iterations
#' @param B2 # of middle loop iterations
#' @param B3 # of innermost loop iterations
#' @param verbose logical variable for printing progress bar. Default to FALSE.
#' @return either 0 or 1 that indicates whether the null hypothesis has been rejected or not
#' @export
single_monitor <- function(y_surveill, yh_bar, sigma, nyears, sig_level, C, B1, B2, B3, verbose=FALSE) {
	fout <- .Call(`singleMonitor`,
		  PACKAGE="TSEnergy",
		  as.vector(y_surveill),
		  as.vector(yh_bar),
		  as.numeric(sigma),
		  as.integer(nyears),
		  as.numeric(sig_level),
		  as.numeric(C),
		  as.integer(B1),
		  as.integer(B2),
		  as.integer(B3),
		  as.logical(verbose))
	fout
}
