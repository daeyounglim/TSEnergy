#include <iostream>
#include <cmath>
#include <RcppArmadillo.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <algorithm>
#include <iterator>
#include <progress.hpp>
#include <progress_bar.hpp>
#include "ListBuilder.h"
#include "random.h"
// [[Rcpp::depends(RcppArmadillo,RcppProgress))]]

// [[Rcpp::export]]
Rcpp::List calc_modelfit_residuals(const arma::vec& yobs,
								   const arma::mat& xobs,
							       const arma::mat& zobs,
							       const arma::mat& uobs,
							       const arma::mat& beta,
							       const arma::mat& alpha,
							       const arma::mat& eta,
							       const arma::mat& sig2,
							       const arma::mat& rho,
							       const arma::mat& lam,
							       const int& nkeep,
							       const bool verbose) {
	double tmp = 0.0;

	mat wobs(T, nkeep, fill::zeros);
	mat xb = xobs * beta;
	{
		Progress prog(nkeep, verbose);
		for (int ikeep = 0; ikeep < nkeep; ++ikeep) {
			vec lam_ikeep = lam.col(ikeep);
			for (int t = 0; t < T; ++t) {
				if (lam_ikeep(t) == 0.0) {
					wobs(t,ikeep) = std::log(yobs(t) - a);
				} else {
					wobs(t,ikeep) = (std::pow(yobs(t) - a, lam_ikeep(t)) - 1.0) / lam_ikeep(t);
				}
			}
		}
		prog.increment();
	}

	vec bayes_resid = arma::mean(wobs - xb, 1);

	return ListBuilder()
	.add("bayes_residuals", bayes_resid);
}
