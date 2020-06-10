#include <iostream>
#include <cmath>
#include <RcppArmadillo.h>
#include <Rmath.h>
#include <Rdefines.h>
#include "quantile.h"
// [[Rcpp::depends(RcppArmadillo))]]

double quantile_cpp(const arma::vec& x, const double& prob) {
	const int n = x.n_elem;
	double idx = 1.0 + static_cast<double>(n-1) * prob;
	int lo = static_cast<int>(std::floor(idx));

	// arma::vec x_sorted = arma::sort(x);
	arma::vec x_sorted = x;
	std::sort(x_sorted.begin(), x_sorted.end());
	double qs = x_sorted(lo-1);
	double g = idx - static_cast<double>(lo);
	return (1.0 - g) * qs + g * x_sorted(lo);
}