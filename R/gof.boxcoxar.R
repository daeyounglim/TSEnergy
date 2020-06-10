"gof.boxcoxar" <- function(object, verbose=FALSE) {
	yobs <- object$yobs
	miss <- object$miss
	a <- object$a
	xobs <- object$xobs
	zobs <- object$zobs
	uobs <- object$uobs
	nkeep <- object$mcmc$nkeep

	pearson.residual <- .Call(`calc_modelfit`,
							 PACKAGE="TSEnergy",
							 as.double(y),
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
	class(gof) <- "gofnmr"
	gof
}