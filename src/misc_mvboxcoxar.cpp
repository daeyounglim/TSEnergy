#include <iostream>
#include <cmath>
#include <RcppArmadillo.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <algorithm>
#include <iterator>
#include "misc_boxcoxar.h"
// [[Rcpp::depends(RcppArmadillo))]]


double loglik_A_prop(const double& zij,
				     const int& index1,
				     const int& index2,
				     const arma::mat& A,
				     const arma::mat& w,
				     const arma::mat& Bx, // = mu
				     const arma::mat& Siginv
				    ) 
{
	const int T = w.n_cols;

	arma::mat A_prop = A;
	A_prop(index1, index2) = (std::exp(2.0 * zij) - 1.0) / (std::exp(2.0 * zij) + 1.0);

	arma::mat resid = w - Bx;
	arma::vec xx = resid.col(0);

	double loglik = -0.5 * arma::dot(xx, Siginv * xx)
					+ 2.0 * zij - 2.0 * arma::log1p(std::exp(2.0 * zij)); // Jacobian for Fisher's z-transformation
	for (int t = 1; t < T; ++t) {
		arma::vec xxx = resid.col(t) - A_prop * w.col(t-1);
		loglik -= 0.5 * arma::dot(xxx, Siginv * xxx);
	}

	return loglik;
}

double loglik_A(const arma::mat& A,
				const arma::mat& w,
				const arma::mat& Bx, // = mu
				const arma::mat& Siginv) 
{
	const int T = w.n_cols;

	arma::mat resid = w - Bx;
	arma::vec xx = resid.col(0);

	double loglik = -0.5 * arma::dot(xx, Siginv * xx)
					+ 2.0 * zij - 2.0 * arma::log1p(std::exp(2.0 * zij)); // Jacobian for Fisher's z-transformation
	for (int t = 1; t < T; ++t) {
		arma::vec xxx = resid.col(t) - A * w.col(t-1);
		loglik -= 0.5 * arma::dot(xxx, Siginv * xxx);
	}

	return loglik;
}

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
				    )
{
	const int K = ystar.n_rows;
	const int T = ystar.n_cols;

	arma::mat C_prop = C;
	C_prop(index1, index2) = cij;
	arma::mat lam = C_prop * xobs.t();

	mat w(arma::size(ystar), arma::fill::zeros); // K * T matrix
    for (int t = 0; t < T; ++t) {
        for (int k = 0; k < K; ++k) {
            if (lam(k,t) == 0.0) {
                w(k,t) = std::log(ystar(k,t) - a);
            } else {
                w(k,t) = (std::pow(ystar(k,t) - a, lam(k,t)) - 1.0) / lam(k,t);
            }
        }
    }

    arma::mat resid = w - Bx;

	arma::vec xx = resid.col(0);
	double loglik = -0.5 * arma::dot(xx, Siginv * xx) + arma::dot(lam - 1.0, arma::log(yobs - a));
	for (int t = 1; t < T; ++t) {
		arma::vec xxx = resid.col(t) - A_prop * w.col(t-1);
		loglik -= 0.5 * arma::dot(xxx, Siginv * xxx);
	}
	return loglik;
}

double loglik_C(const arma::mat& lam,
			    const arma::mat& resid,
			    const arma::mat& w,
			    const arma::mat& A,
			    const arma::mat& Siginv,
			    const double& log_jacobian // log of Jacobian term of the Box-Cox transformation
			    ) {
	const int T = resid.n_cols;

	arma::vec xx = resid.col(0);
	double loglik = -0.5 * arma::dot(xx, Siginv * xx) + log_jacobian;
	for (int t = 1; t < T; ++t) {
		arma::vec xxx = resid.col(t) - A * w.col(t-1);
		loglik -= 0.5 * arma::dot(xxx, Siginv * xxx);
	}
	return loglik;
}
