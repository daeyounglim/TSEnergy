#' get goodness of fit 
#' @param object the output model from fitting a meta analysis/regression model
#' @param type the type of goodness of fit to compute; DIC or LPML or Pearson's residual
#' @param verbose FALSE by default; If TRUE, then progress bar will appear
#' @method gof boxcoxar
#' @export
"gof.boxcoxar" <- function(object, type='pearson', verbose=FALSE) {
	yobs <- object$yobs
	miss <- object$miss
	a <- object$a
	xobs <- object$xobs
	zobs <- object$zobs
	uobs <- object$uobs
	nkeep <- object$mcmc$nkeep
	boxcox_flag <- object$boxcox_flag
	if (boxcox_flag) {
		bayesian.residual <- .Call(`_TSEnergy_calc_modelfit_residuals`,
								 as.double(yobs),
								 as.matrix(xobs),
								 as.matrix(zobs),
								 as.matrix(uobs),
								 as.matrix(object$mcmc.draws$beta),
								 as.matrix(object$mcmc.draws$alpha),
								 as.matrix(object$mcmc.draws$eta),
								 as.matrix(object$mcmc.draws$sig2),
								 as.matrix(object$mcmc.draws$rho),
								 as.matrix(object$mcmc.draws$lam),
								 as.integer(nkeep),
								 as.logical(verbose))
	} else {
		bayesian.residual <- yobs - rowMeans(xobs %*% object$mcmc.draws$beta)
	}
	gof <- list(bayesian.residual = bayesian.residual)
	class(gof) <- "boxcoxar.gof"
	gof
}