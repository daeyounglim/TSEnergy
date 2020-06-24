#ifndef TSENERGY_MISCMVBOXCOXAR_H
#define TSENERGY_MISCMVBOXCOXAR_H
#include <cmath>
#include <Rmath.h>
#include <RcppArmadillo.h>
#include <Rdefines.h>

double loglik_A_prop(const double& zij,
				     const int& index1,
				     const int& index2,
				     const arma::mat& A,
				     const arma::mat& w,
				     const arma::mat& Bx, // = mu
				     const arma::mat& Siginv
				    );
double loglik_A(const arma::mat& A,
				const arma::mat& w,
				const arma::mat& Bx, // = mu
				const arma::mat& Siginv);
double loglik_C_prop(const double& cij,
				     const int& index1,
				     const int& index2,
				     const arma::mat& C,
				     const arma::mat& yobs,
				     const arma::mat& ystar,
				     const arma::mat& xobs,
				     const arma::mat& Bx,
				     const arma::mat& Siginv
				     const double& a
				    );
double loglik_C(const arma::mat& lam,
			    const arma::mat& resid,
			    const arma::mat& w,
			    const arma::mat& A,
			    const arma::mat& Siginv,
			    const double& log_jacobian
			    );

#endif