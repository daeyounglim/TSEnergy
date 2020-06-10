#' get goodness of fit: DIC / LPML / Pearson's residual
#' @param object the output model from fitting a model
#' @param type the type of goodness of fit to compute; DIC / LPML / Pearson's residual
#' @param verbose FALSE by default; If TRUE, then progress bar will appear
#' @return dataframe containing the goodness of fit measure
#' @export
"gof" <- function(object, type='pearson', verbose=FALSE) {
	UseMethod("gof", object)
} 