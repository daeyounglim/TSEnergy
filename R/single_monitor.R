#' @author Daeyoung Lim, \email{daeyoung.lim@uconn.edu}
#' @export
single_monitor <- function(y_surveill, yh_bar, sigma, nyears, sig_level, C, B1, B2, B3, verbose=FALSE) {
	fout <- .Call(`single_monitor`,
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
