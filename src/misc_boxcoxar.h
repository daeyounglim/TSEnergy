#ifndef TSENERGY_MISCBOXCOXAR_H
#define TSENERGY_MISCBOXCOXAR_H
#include <cmath>
#include <Rmath.h>
#include <RcppArmadillo.h>
#include <Rdefines.h>
// [[Rcpp::depends(RcppArmadillo)]]

double loglik_beta(const arma::vec& beta,
			       const arma::vec& wobs,
			       const arma::mat& xobs,
			       const arma::vec& sig2,
			       const double& rho,
			       const double& a,
			       const arma::mat& Siginv,
			       const arma::vec& Siginv_mu);

double loglik_alpha(const arma::vec& alpha,
			 const arma::vec& ystar,
			 const arma::mat& zobs,
			 const arma::vec& sig2,
			 const arma::vec& xb,
			 const double& rho,
			 const double& a,
			 const arma::mat& Siginv,
			 const arma::vec& Siginv_mu);

double loglik_eta(const arma::vec& eta,
				  const arma::vec& wobs,
				  const arma::mat& uobs,
				  const arma::vec& xb,
				  const double& rho,
				  const arma::mat& Siginv,
				  const arma::vec& Siginv_mu);

double loglik_rho(const double& xi,
				  const arma::vec& wobs,
				  const arma::vec& xb,
				  const arma::vec& sig2);

#endif
