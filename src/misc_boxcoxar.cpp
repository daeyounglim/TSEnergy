#include <iostream>
#include <cmath>
#include <RcppArmadillo.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <algorithm>
#include <iterator>
#include "misc_boxcoxar.hpp"
// [[Rcpp::depends(RcppArmadillo))]]
double loglik_beta(const arma::vec& beta,
			       const arma::vec& wobs,
			       const arma::mat& xobs,
			       const arma::vec& sig2,
			       const double& rho,
			       const double& a,
			       const arma::mat& Siginv,
			       const arma::vec& Siginv_mu) {
	using namespace arma;
	vec xb = xobs * beta;
	const int T = wobs.n_elem;
	vec s(T);
	s(0) = wobs(0) - xb(0);
	s.tail(T-1) = wobs.tail(T-1) - xb.tail(T-1) - rho * wobs.head(T-1);
	return -0.5 * accu(arma::square(s) / sig2) -0.5 * arma::dot(beta, Siginv * beta) + arma::dot(beta, Siginv_mu);;
}


double loglik_alpha(const arma::vec& alpha,
					const arma::vec& ystar,
					const arma::mat& zobs,
					const arma::vec& sig2,
					const arma::vec& xb,
					const double& rho,
					const double& a,
					const arma::mat& Siginv,
					const arma::vec& Siginv_mu) {
	using namespace arma;
	vec lam = zobs * alpha;
	const int T = ystar.n_elem;
	vec wobs(T);

	for (int t = 0; t < T; ++t) {
		if (lam(t) == 0.0) {
			wobs(t) = std::log(ystar(t) - a);
		} else {
			wobs(t) = (std::pow(ystar(t) - a, lam(t)) - 1.0) / lam(t);
		}
	}
	vec s(T);
	s(0) = wobs(0) - xb(0);
	s.tail(T-1) = wobs.tail(T-1) - xb.tail(T-1) - rho * wobs.head(T-1);

	return -0.5 * arma::accu(arma::square(s) / sig2) + arma::dot(lam - 1.0, arma::log(ystar - a))
			-0.5 * arma::dot(alpha, Siginv * alpha) + arma::dot(alpha, Siginv_mu);
}


double loglik_eta(const arma::vec& eta,
				  const arma::vec& wobs,
				  const arma::mat& uobs,
				  const arma::vec& xb,
				  const double& rho,
				  const arma::mat& Siginv,
				  const arma::vec& Siginv_mu) {
	using namespace arma;
	vec sig2 = arma::exp(uobs * eta);
	const int T = wobs.n_elem;
	vec s(T);
	s(0) = wobs(0) - xb(0);
	s.tail(T-1) = wobs.tail(T-1) - xb.tail(T-1) - rho * wobs.head(T-1);

	return -0.5 * arma::accu(arma::square(s) / sig2) - 0.5 * accu(arma::log(sig2))
			-0.5 * arma::dot(eta, Siginv * eta) + arma::dot(eta, Siginv_mu);
}

double loglik_rho(const double& xi,
				  const arma::vec& wobs,
				  const arma::vec& xb,
				  const arma::vec& sig2) {
	using namespace arma;
	double rho = (std::exp(2.0 * xi) - 1.0) / (std::exp(2.0 * xi) + 1.0);
	const int T = wobs.n_elem;
	vec s(T-1);
	s = wobs.tail(T-1) - xb.tail(T-1) - rho * wobs.head(T-1);

	return -0.5 * arma::accu(arma::square(s) / sig2.tail(T-1)) + 2.0 * xi - 2.0 * std::log(std::exp(2.0 * xi) + 1.0);
}

