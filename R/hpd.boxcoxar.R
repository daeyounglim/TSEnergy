#' get the highest posterior density (HPD) interval
#' @param object the output model from fitting a meta analysis/regression model
#' @param prob the probability which the HPD interval will cover
#' @return dataframe containing HPD intervals for the parameters

"hpd.boxcoxar" <- function(object, prob = 0.95) {
	out <- list()
	out$beta.hpd <- coda::HPDinterval(mcmc(t(object$mcmc.draws$beta), end=object$mcmc$nkeep), prob=prob)
	out$alpha.hpd <- coda::HPDinterval(mcmc(t(object$mcmc.draws$alpha), end=object$mcmc$nkeep), prob=prob)
	out$eta.hpd <- coda::HPDinterval(mcmc(t(object$mcmc.draws$eta), end=object$mcmc$nkeep), prob=prob)
	out$sig2.hpd <- coda::HPDinterval(mcmc(t(object$mcmc.draws$sig2), end=object$mcmc$nkeep), prob=prob)
	out$lam.hpd <- coda::HPDinterval(mcmc(t(object$mcmc.draws$lam), end=object$mcmc$nkeep), prob=prob)
	
	class(out) <- "boxcoxar.hpd"
	out
}