#' get goodness of fit 
#' @param object the output model from fitting a model
#' @param type the type of goodness of fit to compute; DIC / LPML / Pearson's residual
#' @param verbose FALSE by default; If TRUE, then progress bar will appear

"gof.boxcoxar" <- function(object, type='pearson', verbose=FALSE) {
	yobs <- object$yobs
	miss <- object$miss
	a <- object$a
	xobs <- object$xobs
	zobs <- object$zobs
	uobs <- object$uobs
	nkeep <- object$mcmc$nkeep

	pearson.residual <- .Call(`calc_modelfit`,
							 PACKAGE="TSEnergy",
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
	class(gof) <- "boxcoxar.gof"
	gof
}