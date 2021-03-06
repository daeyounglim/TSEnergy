#' get the highest posterior density (HPD) interval
#' @param object the output model from fitting a model
#' @param prob the probability which the HPD interval will cover
#' @return dataframe containing HPD intervals for the parameters
#' @export
"hpd" <- function(object, prob) {
	UseMethod("hpd", object)
} 