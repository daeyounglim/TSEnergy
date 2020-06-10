#include <iostream>
#include <cmath>
#include <RcppArmadillo.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <algorithm>
#include <iterator>
#include <progress.hpp>
#include <progress_bar.hpp>
#include "quantile.h"
#include "ListBuilder.h"
// [[Rcpp::depends(RcppArmadillo,RcppProgress))]]

int single_monitor(const arma::vec& y_surveill,
						 const arma::vec& yh_bar, // monthly average of historical data
						 const double& sigma,
						 const double& nyears, // number of historical years (= "T" in the paper)
						 const double& sig_level, // significance level
						 const double& c, // alternative H_1 : Y_k ~ N(c * beta, sig2)
						 const int& B1,
						 const int& B2,
						 const int& B3,
						 const bool verbose) {
	using namespace arma;
	using namespace Rcpp;
	using namespace R;
	using namespace std;

	const int nmonths = 12;
	double z_alpha = qnorm5(1.0-sig_level, 0.0, 1.0, 1, 0);
	vec pvalue_hat(B2, fill::zeros);
	for (int b2=0; b2 < B2; ++b2) {
		/* generate fake data */
		vec y_tilde(nmonths);
		std::generate(y_tilde.begin(), y_tilde.end(), ::norm_rand);
		y_tilde *= sigma;
		y_tilde += yh_bar;

		vec delta(B1, fill::zeros);
		for (int b1=0; b1 < B1; ++b1) {
			/* sample from fiducial distribution */
			vec beta_fiduc(nmonths);
			std::generate(beta_fiduc.begin(), beta_fiduc.end(), ::norm_rand);
			beta_fiduc *= sigma / std::sqrt(nyears);
			beta_fiduc += yh_bar;

			double sig_fiduc = std::sqrt(6.0 * (nyears - 1.0) * std::pow(sigma, 2.0) / ::Rf_rgamma(6.0 * (nyears-1.0) ,1.0));
			double std_test = arma::accu(beta_fiduc % (y_tilde - beta_fiduc) / (std::sqrt(sig_fiduc) * arma::norm(beta_fiduc)));
			delta(b1) = (std_test > z_alpha) ? 1.0 : 0.0;
		}
		pvalue_hat(b2) = arma::mean(delta);
	}
	double gamma_alpha = quantile_cpp(pvalue_hat, sig_level);
	vec delta(B3, fill::zeros);
	for (int b=0; b < B3; ++b) {
		/* sample from fiducial distribution */
		vec beta_fiduc(nmonths);
		std::generate(beta_fiduc.begin(), beta_fiduc.end(), ::norm_rand);
		beta_fiduc *= sigma / std::sqrt(nyears);
		beta_fiduc += yh_bar;

		double sig_fiduc = std::sqrt(6.0 * (nyears - 1.0) * std::pow(sigma, 2.0) / ::Rf_rgamma(6.0 * (nyears-1.0) ,1.0));
		double std_test = arma::accu(beta_fiduc % (y_surveill - beta_fiduc) / (std::sqrt(sig_fiduc) * arma::norm(beta_fiduc)));
		delta(b) = (std_test > z_alpha) ? 1.0 : 0.0;
	}
	double p_surveill = arma::mean(delta);
	return (p_surveill >= gamma_alpha) ? 1 : 0;
}
