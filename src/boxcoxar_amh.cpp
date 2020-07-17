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
#include "misc_boxcoxar.h"
#include "nelmin.h"
// [[Rcpp::depends(RcppArmadillo,RcppProgress))]]

// [[Rcpp::export]]
Rcpp::List boxcoxar_amh(const arma::vec& yobs,
					const arma::uvec& miss,
					const double& a,
				    const arma::mat& xobs,
				    const arma::mat& zobs,
				    const arma::mat& uobs,
				    const arma::vec& mu_beta0,
				    const arma::mat& Sig_beta0,
				    const arma::vec& mu_alpha0,
				    const arma::mat& Sig_alpha0,
				    const arma::vec& mu_eta0,
				    const arma::mat& Sig_eta0,
				    const int& ndiscard,
				    const int& nskip,
				    const int& nkeep,
				    const bool boxcox_flag,
				    const bool verbose) {
	using namespace arma;
	using namespace Rcpp;
	using namespace R;
	using namespace std;

	const int T = yobs.n_elem;
	const int xcols = xobs.n_cols;
	const int zcols = zobs.n_cols;
	const int ucols = uobs.n_cols;

	// initialize parameters
	vec ystar = yobs;
	vec beta(xcols, fill::zeros);
	vec alpha(zcols);
	alpha.fill(0.01);
	vec eta(ucols, fill::zeros);
	vec lam = zobs * alpha;
	double rho = 0.5;
	vec sig2 = arma::exp(uobs * eta);
	vec wobs(T, fill::zeros);
	if (boxcox_flag) {
		for (int t = 0; t < T; ++t) {
			if (lam(t) == 0.0) {
				wobs(t) = std::log(ystar(t) - a);
			} else {
				wobs(t) = (std::pow(ystar(t) - a, lam(t)) - 1.0) / lam(t);
			}
		}
	} else {
		wobs = ystar;
	}
	bool miss_flag = any(miss);
	uvec miss_idx = find(miss);
	if (miss_flag) {
		for (unsigned int t = 0; t < miss_idx.n_elem; ++t) {
			int idx = miss_idx(t);
			double normrv = ::norm_rand();
			ystar(idx) = normrv;
			if (boxcox_flag) {
				if (lam(idx) == 0.0) {
					wobs(idx) = std::log(ystar(idx) - a);
				} else {
					wobs(idx) = (std::pow(ystar(idx) - a, lam(idx)) - 1.0) / lam(idx);
				}
			} else {
				wobs(idx) = ystar(idx);
			}
		}
	}

	// Miscellaneous
	vec xb = xobs * beta;
	mat Siginv_beta0 = inv(Sig_beta0);
	mat Siginv_alpha0 = inv(Sig_alpha0);
	mat Siginv_eta0 = inv(Sig_eta0);

	vec Siginv_mu_beta = Siginv_beta0 * mu_beta0;
	vec Siginv_mu_alpha = Siginv_alpha0 * mu_alpha0;
	vec Siginv_mu_eta = Siginv_eta0 * mu_eta0;

	/*************************
	Parameters for adaptive MH
	*************************/
	const double obj_rate = 0.44;
	const int batch_length = std::min(50, ndiscard+nkeep*nskip);
	const int batch_total = (ndiscard + nskip*nkeep) / batch_length;

	vec alpha_accepts(zcols, fill::zeros);
	mat alpha_rates(zcols, batch_total, fill::zeros);
	vec alpha_tunings(zcols);
	alpha_tunings.fill(0.1);
	int batch_num = 0;

	vec eta_accepts(ucols, fill::zeros);
	mat eta_rates(ucols, batch_total, fill::zeros);
	vec eta_tunings(ucols);
	eta_tunings.fill(0.1);

	vec beta_accepts(xcols, fill::zeros);
	mat beta_rates(xcols, batch_total, fill::zeros);
	vec beta_tunings(xcols);
	beta_tunings.fill(0.1);

	double rho_accepts = 0.0;
	vec rho_rates(batch_total, fill::zeros);
	double rho_tunings = 0.1;


	/*******************
	Begin burn-in period
	*******************/
	if (verbose) {
		Rcout << "Warming up" << endl;
	}
	{
		Progress prog(ndiscard, verbose);
		for (int idiscard = 0; idiscard < ndiscard; ++idiscard) {
			if (Progress::check_abort()) {
				return Rcpp::List::create(Rcpp::Named("error") = "user interrupt aborted");
			}
			/**********
			Sample beta
			**********/
			for (int j = 0; j < xcols; ++j) {
				vec beta_prop = beta;
				beta_prop(j) = ::norm_rand() * std::exp(beta_tunings(j)) + beta(j);
				double ll_diff = loglik_beta(beta_prop, wobs, xobs, sig2, rho, a, Siginv_beta0, Siginv_mu_beta)
								- loglik_beta(beta, wobs, xobs, sig2, rho, a, Siginv_beta0, Siginv_mu_beta);
				if (std::log(::unif_rand()) < ll_diff) {
					beta(j) = beta_prop(j);
					++beta_accepts(j);
				}
			}
			xb = xobs * beta;

			/***********
			Sample alpha
			***********/
			if (boxcox_flag) {
				for (int j = 0; j < zcols; ++j) {
					vec alpha_prop = alpha;
					alpha_prop(j) = ::norm_rand() * std::exp(alpha_tunings(j)) + alpha(j);

					// log-likelihood difference
					double ll_diff = loglik_alpha(alpha_prop, ystar, zobs, sig2, xb, rho, a, Siginv_alpha0, Siginv_mu_alpha) -
							  loglik_alpha(alpha, ystar, zobs, sig2, xb, rho, a, Siginv_alpha0, Siginv_mu_alpha);
					if (std::log(::unif_rand()) < ll_diff) {
						alpha(j) = alpha_prop(j);
						++alpha_accepts(j);
					}
				}
				lam = zobs * alpha;
				for (int t = 0; t < T; ++t) {
					if (lam(t) == 0.0) {
						wobs(t) = std::log(ystar(t) - a);
					} else {
						wobs(t) = (std::pow(ystar(t) - a, lam(t)) - 1.0) / lam(t);
					}
				}
			}


			/*********
			Sample eta
			*********/
			for (int j = 0; j < ucols; ++j) {
				vec eta_prop = eta;
				eta_prop(j) = ::norm_rand() * std::exp(eta_tunings(j)) + eta(j);

				// log-likelihood difference
				double ll_diff = loglik_eta(eta_prop, wobs, uobs, xb, rho, Siginv_eta0, Siginv_mu_eta) -
						  loglik_eta(eta, wobs, uobs, xb, rho, Siginv_eta0, Siginv_mu_eta);
				if (std::log(::unif_rand()) < ll_diff) {
					eta(j) = eta_prop(j);
					++eta_accepts(j);
				}
			}
			sig2 = arma::exp(uobs * eta);

			/*********
			Sample rho
			*********/
			double xi = 0.5 * std::log((1.0 + rho) / (1.0 - rho));
			double xi_prop = ::norm_rand() * std::exp(rho_tunings) + xi;

			// log-likelihood difference
			double ll_diff = loglik_rho(xi_prop, wobs, xb, sig2) -
							 loglik_rho(xi, wobs, xb, sig2);
			if (std::log(::unif_rand()) < ll_diff) {
				rho = (std::exp(2.0 * xi_prop) - 1.0) / (std::exp(2.0 * xi_prop) + 1.0);
				++rho_accepts;
			}

			if ((idiscard+1) % batch_length == 0) {
				rho_accepts /= static_cast<double>(batch_length);
				rho_rates(batch_num) = rho_accepts;
				if (rho_accepts > obj_rate) {
					rho_tunings += std::min(0.01, 1.0 / std::sqrt(static_cast<double>(batch_num+1)));
				} else {
					rho_tunings -= std::min(0.01, 1.0 / std::sqrt(static_cast<double>(batch_num+1)));
				}

				beta_accepts /= static_cast<double>(batch_length);
				beta_rates.col(batch_num) = beta_accepts;
				for (int i=0; i<xcols; ++i) {
					if (beta_accepts(i) > obj_rate) {
						beta_tunings(i) += std::min(0.01, 1.0 / std::sqrt(static_cast<double>(batch_num+1)));
					} else {
						beta_tunings(i) -= std::min(0.01, 1.0 / std::sqrt(static_cast<double>(batch_num+1)));
					}
				}

				if (boxcox_flag) {
					alpha_accepts /= static_cast<double>(batch_length);
					alpha_rates.col(batch_num) = alpha_accepts;
					for (int i=0; i<zcols; ++i) {
						if (alpha_accepts(i) > obj_rate) {
							alpha_tunings(i) += std::min(0.01, 1.0 / std::sqrt(static_cast<double>(batch_num+1)));
						} else {
							alpha_tunings(i) -= std::min(0.01, 1.0 / std::sqrt(static_cast<double>(batch_num+1)));
						}
					}
				}

				eta_accepts /= static_cast<double>(batch_length);
				eta_rates.col(batch_num) = eta_accepts;
				for (int i=0; i<ucols; ++i) {
					if (eta_accepts(i) > obj_rate) {
						eta_tunings(i) += std::min(0.01, 1.0 / std::sqrt(static_cast<double>(batch_num+1)));
					} else {
						eta_tunings(i) -= std::min(0.01, 1.0 / std::sqrt(static_cast<double>(batch_num+1)));
					}
				}
				++batch_num;
				rho_accepts = 0.0;
				beta_accepts.fill(0.0);
				if (boxcox_flag) {
					alpha_accepts.fill(0.0);
				}
				eta_accepts.fill(0.0);
			}
			prog.increment();
		}
	}

	mat beta_save(xcols, nkeep, fill::zeros);
	mat alpha_save(zcols, nkeep, fill::zeros);
	mat lam_save(T, nkeep, fill::zeros);
	mat eta_save(ucols, nkeep, fill::zeros);
	mat sig2_save(T, nkeep, fill::zeros);
	vec rho_save(nkeep, fill::zeros);
	/*******************
	Begin sampling period
	*******************/
	int icount_mh = ndiscard;
	if (verbose) {
		Rcout << "Sampling posterior" << endl;
	}
	{
		Progress prog(nkeep, verbose);
		for (int ikeep = 0; ikeep < nkeep; ++ikeep) {
			if (Progress::check_abort()) {
				return Rcpp::List::create(Rcpp::Named("error") = "user interrupt aborted");
			}
			for (int iskip = 0; iskip < nskip; ++iskip) {
				++icount_mh;
				/**********
				Sample beta
				**********/
				for (int j = 0; j < xcols; ++j) {
					vec beta_prop = beta;
					beta_prop(j) = ::norm_rand() * std::exp(beta_tunings(j)) + beta(j);
					double ll_diff = loglik_beta(beta_prop, wobs, xobs, sig2, rho, a, Siginv_beta0, Siginv_mu_beta)
									- loglik_beta(beta, wobs, xobs, sig2, rho, a, Siginv_beta0, Siginv_mu_beta);
					if (std::log(::unif_rand()) < ll_diff) {
						beta(j) = beta_prop(j);
						++beta_accepts(j);
					}
				}
				xb = xobs * beta;

				if (boxcox_flag) {
					/***********
					Sample alpha
					***********/
					for (int j = 0; j < zcols; ++j) {
						vec alpha_prop = alpha;
						alpha_prop(j) = ::norm_rand() * std::exp(alpha_tunings(j)) + alpha(j);

						// log-likelihood difference
						double ll_diff = loglik_alpha(alpha_prop, ystar, zobs, sig2, xb, rho, a, Siginv_alpha0, Siginv_mu_alpha) -
								  loglik_alpha(alpha, ystar, zobs, sig2, xb, rho, a, Siginv_alpha0, Siginv_mu_alpha);
						if (std::log(::unif_rand()) < ll_diff) {
							alpha(j) = alpha_prop(j);
							++alpha_accepts(j);
						}
					}
					lam = zobs * alpha;
					for (int t = 0; t < T; ++t) {
						if (lam(t) == 0.0) {
							wobs(t) = std::log(ystar(t) - a);
						} else {
							wobs(t) = (std::pow(ystar(t) - a, lam(t)) - 1.0) / lam(t);
						}
					}
				}

				/*********
				Sample eta
				*********/
				for (int j = 0; j < ucols; ++j) {
					vec eta_prop = eta;
					eta_prop(j) = ::norm_rand() * std::exp(eta_tunings(j)) + eta(j);

					// log-likelihood difference
					double ll_diff = loglik_eta(eta_prop, wobs, uobs, xb, rho, Siginv_eta0, Siginv_mu_eta) -
							  loglik_eta(eta, wobs, uobs, xb, rho, Siginv_eta0, Siginv_mu_eta);
					if (std::log(::unif_rand()) < ll_diff) {
						eta(j) = eta_prop(j);
						++eta_accepts(j);
					}
				}
				sig2 = arma::exp(uobs * eta);

				/*********
				Sample rho
				*********/
				double xi = 0.5 * std::log((1.0 + rho) / (1.0 - rho));
				double xi_prop = ::norm_rand() * std::exp(rho_tunings) + xi;

				// log-likelihood difference
				double ll_diff = loglik_rho(xi_prop, wobs, xb, sig2) -
								 loglik_rho(xi, wobs, xb, sig2);
				if (std::log(::unif_rand()) < ll_diff) {
					rho = (std::exp(2.0 * xi_prop) - 1.0) / (std::exp(2.0 * xi_prop) + 1.0);
					++rho_accepts;
				}

				if (icount_mh % batch_length == 0) {
					rho_accepts /= static_cast<double>(batch_length);
					rho_rates(batch_num) = rho_accepts;
					if (rho_accepts > obj_rate) {
						rho_tunings += std::min(0.01, 1.0 / std::sqrt(static_cast<double>(batch_num+1)));
					} else {
						rho_tunings -= std::min(0.01, 1.0 / std::sqrt(static_cast<double>(batch_num+1)));
					}

					beta_accepts /= static_cast<double>(batch_length);
					beta_rates.col(batch_num) = beta_accepts;
					for (int i=0; i<xcols; ++i) {
						if (beta_accepts(i) > obj_rate) {
							beta_tunings(i) += std::min(0.01, 1.0 / std::sqrt(static_cast<double>(batch_num+1)));
						} else {
							beta_tunings(i) -= std::min(0.01, 1.0 / std::sqrt(static_cast<double>(batch_num+1)));
						}
					}

					if (boxcox_flag) {
						alpha_accepts /= static_cast<double>(batch_length);
						alpha_rates.col(batch_num) = alpha_accepts;
						for (int i=0; i<zcols; ++i) {
							if (alpha_accepts(i) > obj_rate) {
								alpha_tunings(i) += std::min(0.01, 1.0 / std::sqrt(static_cast<double>(batch_num+1)));
							} else {
								alpha_tunings(i) -= std::min(0.01, 1.0 / std::sqrt(static_cast<double>(batch_num+1)));
							}
						}
					}

					eta_accepts /= static_cast<double>(batch_length);
					eta_rates.col(batch_num) = eta_accepts;
					for (int i=0; i<ucols; ++i) {
						if (eta_accepts(i) > obj_rate) {
							eta_tunings(i) += std::min(0.01, 1.0 / std::sqrt(static_cast<double>(batch_num+1)));
						} else {
							eta_tunings(i) -= std::min(0.01, 1.0 / std::sqrt(static_cast<double>(batch_num+1)));
						}
					}
					++batch_num;
					rho_accepts = 0.0;
					if (boxcox_flag) {
						alpha_accepts.fill(0.0);
					}
					beta_accepts.fill(0.0);
					eta_accepts.fill(0.0);
				}
			}
			beta_save.col(ikeep) = beta;
			if (boxcox_flag) {
				alpha_save.col(ikeep) = alpha;
			}
			eta_save.col(ikeep) = eta;
			lam_save.col(ikeep) = lam;
			sig2_save.col(ikeep) = sig2;
			rho_save(ikeep) = rho;
			prog.increment();
		}
	}
	
	if (boxcox_flag) {
		return ListBuilder()
		.add("beta", beta_save)
		.add("alpha", alpha_save)
		.add("eta", eta_save)
		.add("sig2", sig2_save)
		.add("rho", rho_save)
		.add("lam", lam_save)
		.add("beta_acceptance", beta_rates)
		.add("alpha_acceptance", alpha_rates)
		.add("eta_acceptance", eta_rates)
		.add("rho_acceptance", rho_rates);
	} else {
		return ListBuilder()
		.add("beta", beta_save)
		.add("eta", eta_save)
		.add("sig2", sig2_save)
		.add("rho", rho_save)
		.add("beta_acceptance", beta_rates)
		.add("alpha_acceptance", alpha_rates)
		.add("eta_acceptance", eta_rates)
		.add("rho_acceptance", rho_rates);
	}
}