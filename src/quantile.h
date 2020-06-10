#ifndef TSEnergy_QUANTILE_CPP_H
#define TSEnergy_QUANTILE_CPP_H
#include <iostream>
#include <cmath>
#include <RcppArmadillo.h>
#include <Rmath.h>
#include <Rdefines.h>
// [[Rcpp::depends(RcppArmadillo))]]

double quantile_cpp(const arma::vec& x, const double& prob);

#endif
