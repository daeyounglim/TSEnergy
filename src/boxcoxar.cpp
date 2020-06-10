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

Rcpp::List boxcoxar(const arma::vec& yobs,
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
	for (int t = 0; t < T; ++t) {
		if (lam(t) == 0.0) {
			wobs(t) = std::log(ystar(t) - a);
		} else {
			wobs(t) = (std::pow(ystar(t) - a, lam(t)) - 1.0) / lam(t);
		}
	}
	bool miss_flag = any(miss);
	uvec miss_idx = find(miss);
	if (miss_flag) {
		for (unsigned int t = 0; t < miss_idx.n_elem; ++t) {
			int idx = miss_idx(t);
			double normrv = ::norm_rand();
			ystar(idx) = normrv;
			if (lam(idx) == 0.0) {
				wobs(idx) = std::log(ystar(idx) - a);
			} else {
				wobs(idx) = (std::pow(ystar(idx) - a, lam(idx)) - 1.0) / lam(idx);
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

	vec alpha_rates(zcols, fill::zeros);
	vec eta_rates(ucols, fill::zeros);
	double rho_rates = 0.0;


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
			mat Ainv_beta = xobs.t() * diagmat(1.0 / sig2) * xobs + Siginv_beta0;
			mat Ainv_chol_beta = chol(Ainv_beta);
			vec a_beta_tmp(T);
			a_beta_tmp(0) = wobs(0) / sig2(0);
			a_beta_tmp.tail(T-1) = (wobs.tail(T-1) - rho * wobs.head(T-1)) / sig2.tail(T-1);

			vec a_beta = arma::solve(arma::trimatu(Ainv_chol_beta), arma::solve(trimatl(Ainv_chol_beta.t()), xobs.t() * a_beta_tmp + Siginv_mu_beta));
			vec btmp(xcols);
			std::generate(btmp.begin(), btmp.end(), ::norm_rand);
			beta = a_beta + arma::solve(arma::trimatu(Ainv_chol_beta), btmp);
			xb = xobs * beta;

			/***********
			Sample alpha
			***********/
			for (int j = 0; j < zcols; ++j) {
				vec alpha_prop = alpha;
				auto fx = [&](double alphaj[])->double {
					alpha_prop(j) = alphaj[0];
					return -loglik_alpha(alpha_prop, ystar, zobs, sig2, xb, rho, a, Siginv_alpha0, Siginv_mu_alpha);
				};
				double start[] = { alpha(j) };
				double xmin[] = { 0.0 };
				double ynewlo = 0.0;
				double reqmin = 1.0e-20;
				int konvge = 5;
				int kcount = 1000;
				double step[] = { 0.02 };
				int icount = 0;
				int numres = 0;
				int ifault = 0;
				nelmin(fx, 1, start, xmin, &ynewlo, reqmin, step, konvge, kcount, &icount, &numres, &ifault);

				double xmax = xmin[0];
				double minll = ynewlo;

				mat cl(5,3, fill::zeros);
				vec dl(5, fill::zeros);
				double step_size = 0.2;
				alpha_burnin_block:
					for (int iii=0; iii < 5; ++iii) {
						double e1 = static_cast<double>(iii-2);
						cl(iii,0) = std::pow(xmax + e1 * step_size, 2.0);
						cl(iii,1) = xmax + e1 * step_size;
						cl(iii,2) = 1.0;
						alpha_prop(j) = xmax + e1 * step_size;
						dl(iii) = -loglik_alpha(alpha_prop, ystar, zobs, sig2, xb, rho, a, Siginv_alpha0, Siginv_mu_alpha);
					}

				for (int ni=0; ni < 5; ++ni) {
					if ((ni+1) != 3) {
						if (dl(ni) <= minll) {
							step_size *= 1.2;
							goto alpha_burnin_block;
						}
					}
				}

				vec fl = solve(cl.t() * cl, cl.t() * dl);
				double sigmaa = std::sqrt(0.5 / fl(0));

				alpha_prop(j) = ::norm_rand() * sigmaa + xmax;

				// log-likelihood difference
				double ll_diff = loglik_alpha(alpha_prop, ystar, zobs, sig2, xb, rho, a, Siginv_alpha0, Siginv_mu_alpha) -
						  loglik_alpha(alpha, ystar, zobs, sig2, xb, rho, a, Siginv_alpha0, Siginv_mu_alpha) -
						  0.5 * (std::pow(alpha(j) - xmax, 2.0) - std::pow(alpha_prop(j) - xmax, 2.0)) / std::pow(sigmaa, 2.0);
				if (std::log(::unif_rand()) < ll_diff) {
					alpha(j) = alpha_prop(j);
					++alpha_rates(j);
				}
			}
			lam = zobs * alpha;

			/*********
			Sample eta
			*********/
			for (int j = 0; j < ucols; ++j) {
				vec eta_prop = eta;
				auto fx = [&](double etaj[])->double {
					eta_prop(j) = etaj[0];
					return -loglik_eta(eta_prop, wobs, uobs, xb, rho, Siginv_eta0, Siginv_mu_eta);
				};
				double start[] = { eta(j) };
				double xmin[] = { 0.0 };
				double ynewlo = 0.0;
				double reqmin = 1.0e-20;
				int konvge = 5;
				int kcount = 1000;
				double step[] = { 0.05 };
				int icount = 0;
				int numres = 0;
				int ifault = 0;
				nelmin(fx, 1, start, xmin, &ynewlo, reqmin, step, konvge, kcount, &icount, &numres, &ifault);

				double xmax = xmin[0];
				double minll = ynewlo;

				mat cl(5,3, fill::zeros);
				vec dl(5, fill::zeros);
				double step_size = 0.2;
				eta_burnin_block:
					for (int iii=0; iii < 5; ++iii) {
						double e1 = static_cast<double>(iii-2);
						cl(iii,0) = std::pow(xmax + e1 * step_size, 2.0);
						cl(iii,1) = xmax + e1 * step_size;
						cl(iii,2) = 1.0;
						eta_prop(j) = xmax + e1 * step_size;
						dl(iii) = -loglik_eta(eta_prop, wobs, uobs, xb, rho, Siginv_eta0, Siginv_mu_eta);
					}

				for (int ni=0; ni < 5; ++ni) {
					if ((ni+1) != 3) {
						if (dl(ni) <= minll) {
							step_size *= 1.2;
							goto eta_burnin_block;
						}
					}
				}

				vec fl = solve(cl.t() * cl, cl.t() * dl);
				double sigmaa = std::sqrt(0.5 / fl(0));

				eta_prop(j) = ::norm_rand() * sigmaa + xmax;

				// log-likelihood difference
				double ll_diff = loglik_eta(eta_prop, wobs, uobs, xb, rho, Siginv_eta0, Siginv_mu_eta) -
						  loglik_eta(eta, wobs, uobs, xb, rho, Siginv_eta0, Siginv_mu_eta) -
						  0.5 * (std::pow(eta(j) - xmax, 2.0) - std::pow(eta_prop(j) - xmax, 2.0)) / std::pow(sigmaa, 2.0);
				if (std::log(::unif_rand()) < ll_diff) {
					eta(j) = eta_prop(j);
					++eta_rates(j);
				}
			}
			sig2 = arma::exp(uobs * eta);

			/*********
			Sample rho
			*********/
			double xi = 0.5 * std::log((1.0 + rho) / (1.0 - rho));
			auto fx_rho = [&](double xi_input[])->double {
				return -loglik_rho(xi_input[0], wobs, xb, sig2);
			};
			double start[] = { xi };
			double xmin[] = { 0.0 };
			double ynewlo = 0.0;
			double reqmin = 1.0e-20;
			int konvge = 5;
			int kcount = 1000;
			double step[] = { 0.02 };
			int icount = 0;
			int numres = 0;
			int ifault = 0;
			nelmin(fx_rho, 1, start, xmin, &ynewlo, reqmin, step, konvge, kcount, &icount, &numres, &ifault);
			double xmax = xmin[0];
			double minll = ynewlo;

			mat cl(5,3, fill::zeros);
			vec dl(5, fill::zeros);
			double step_size = 0.2;
			rho_burnin_block:
				for (int iii=0; iii < 5; ++iii) {
					double e1 = static_cast<double>(iii-2);
					cl(iii,0) = std::pow(xmax + e1 * step_size, 2.0);
					cl(iii,1) = xmax + e1 * step_size;
					cl(iii,2) = 1.0;
					dl(iii) = -loglik_rho(xmax + e1 * step_size, wobs, xb, sig2);
				}

			for (int ni=0; ni < 5; ++ni) {
				if ((ni+1) != 3) {
					if (dl(ni) <= minll) {
						step_size *= 1.2;
						goto rho_burnin_block;
					}
				}
			}

			vec fl = arma::solve(cl.t() * cl, cl.t() * dl);
			double sigmaa = std::sqrt(0.5 / fl(0));


			double xi_prop = ::norm_rand() * sigmaa + xmax;

			// log-likelihood difference
			double ll_diff = loglik_rho(xi_prop, wobs, xb, sig2) -
							 loglik_rho(xi, wobs, xb, sig2) -
							 0.5 * (std::pow(xi - xmax, 2.0) - std::pow(xi_prop - xmax, 2.0)) / std::pow(sigmaa, 2.0);
			if (std::log(::unif_rand()) < ll_diff) {
				rho = (std::exp(2.0 * xi_prop) - 1.0) / (std::exp(2.0 * xi_prop) + 1.0);
				++rho_rates;
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
				/**********
				Sample beta
				**********/
				mat Ainv_beta = xobs.t() * diagmat(1.0 / sig2) * xobs + Siginv_beta0;
				mat Ainv_chol_beta = chol(Ainv_beta);
				vec a_beta_tmp(T);
				a_beta_tmp(0) = wobs(0) / sig2(0);
				a_beta_tmp.tail(T-1) = (wobs.tail(T-1) - rho * wobs.head(T-1)) / sig2.tail(T-1);

				vec a_beta = arma::solve(arma::trimatu(Ainv_chol_beta), arma::solve(trimatl(Ainv_chol_beta.t()), xobs.t() * a_beta_tmp + Siginv_mu_beta));
				vec btmp(xcols);
				std::generate(btmp.begin(), btmp.end(), ::norm_rand);
				beta = a_beta + arma::solve(arma::trimatu(Ainv_chol_beta), btmp);
				xb = xobs * beta;

				/***********
				Sample alpha
				***********/
				for (int j = 0; j < zcols; ++j) {
					vec alpha_prop = alpha;
					auto fx = [&](double alphaj[])->double {
						alpha_prop(j) = alphaj[0];
						return -loglik_alpha(alpha_prop, ystar, zobs, sig2, xb, rho, a, Siginv_alpha0, Siginv_mu_alpha);
					};
					double start[] = { alpha(j) };
					double xmin[] = { 0.0 };
					double ynewlo = 0.0;
					double reqmin = 1.0e-20;
					int konvge = 5;
					int kcount = 1000;
					double step[] = { 0.02 };
					int icount = 0;
					int numres = 0;
					int ifault = 0;
					nelmin(fx, 1, start, xmin, &ynewlo, reqmin, step, konvge, kcount, &icount, &numres, &ifault);

					double xmax = xmin[0];
					double minll = ynewlo;

					mat cl(5,3, fill::zeros);
					vec dl(5, fill::zeros);
					double step_size = 0.1;
					alpha_sampling_block:
						for (int iii=0; iii < 5; ++iii) {
							double e1 = static_cast<double>(iii-2);
							cl(iii,0) = std::pow(xmax + e1 * step_size, 2.0);
							cl(iii,1) = xmax + e1 * step_size;
							cl(iii,2) = 1.0;
							alpha_prop(j) = xmax + e1 * step_size;
							dl(iii) = -loglik_alpha(alpha_prop, ystar, zobs, sig2, xb, rho, a, Siginv_alpha0, Siginv_mu_alpha);
						}

					for (int ni=0; ni < 5; ++ni) {
						if ((ni+1) != 3) {
							if (dl(ni) <= minll) {
								step_size *= 1.2;
								goto alpha_sampling_block;
							}
						}
					}

					vec fl = solve(cl.t() * cl, cl.t() * dl);
					double sigmaa = std::sqrt(0.5 / fl(0));

					alpha_prop(j) = ::norm_rand() * sigmaa + xmax;

					// log-likelihood difference
					double ll_diff = loglik_alpha(alpha_prop, ystar, zobs, sig2, xb, rho, a, Siginv_alpha0, Siginv_mu_alpha) -
							  loglik_alpha(alpha, ystar, zobs, sig2, xb, rho, a, Siginv_alpha0, Siginv_mu_alpha) -
							  0.5 * (std::pow(alpha(j) - xmax, 2.0) - std::pow(alpha_prop(j) - xmax, 2.0)) / std::pow(sigmaa, 2.0);
					if (std::log(::unif_rand()) < ll_diff) {
						alpha(j) = alpha_prop(j);
						++alpha_rates(j);
					}
				}
				lam = zobs * alpha;

				/*********
				Sample eta
				*********/
				for (int j = 0; j < ucols; ++j) {
					vec eta_prop = eta;
					auto fx = [&](double etaj[])->double {
						eta_prop(j) = etaj[0];
						return -loglik_eta(eta_prop, wobs, uobs, xb, rho, Siginv_eta0, Siginv_mu_eta);
					};
					double start[] = { eta(j) };
					double xmin[] = { 0.0 };
					double ynewlo = 0.0;
					double reqmin = 1.0e-20;
					int konvge = 5;
					int kcount = 1000;
					double step[] = { 0.05 };
					int icount = 0;
					int numres = 0;
					int ifault = 0;
					nelmin(fx, 1, start, xmin, &ynewlo, reqmin, step, konvge, kcount, &icount, &numres, &ifault);

					double xmax = xmin[0];
					double minll = ynewlo;

					mat cl(5,3, fill::zeros);
					vec dl(5, fill::zeros);
					double step_size = 0.1;
					eta_sampling_block:
						for (int iii=0; iii < 5; ++iii) {
							double e1 = static_cast<double>(iii-2);
							cl(iii,0) = std::pow(xmax + e1 * step_size, 2.0);
							cl(iii,1) = xmax + e1 * step_size;
							cl(iii,2) = 1.0;
							eta_prop(j) = xmax + e1 * step_size;
							dl(iii) = -loglik_eta(eta_prop, wobs, uobs, xb, rho, Siginv_eta0, Siginv_mu_eta);
						}

					for (int ni=0; ni < 5; ++ni) {
						if ((ni+1) != 3) {
							if (dl(ni) <= minll) {
								step_size *= 1.2;
								goto eta_sampling_block;
							}
						}
					}

					vec fl = solve(cl.t() * cl, cl.t() * dl);
					double sigmaa = std::sqrt(0.5 / fl(0));

					eta_prop(j) = ::norm_rand() * sigmaa + xmax;

					// log-likelihood difference
					double ll_diff = loglik_eta(eta_prop, wobs, uobs, xb, rho, Siginv_eta0, Siginv_mu_eta) -
							  loglik_eta(eta, wobs, uobs, xb, rho, Siginv_eta0, Siginv_mu_eta) -
							  0.5 * (std::pow(eta(j) - xmax, 2.0) - std::pow(eta_prop(j) - xmax, 2.0)) / std::pow(sigmaa, 2.0);
					if (std::log(::unif_rand()) < ll_diff) {
						eta(j) = eta_prop(j);
						++eta_rates(j);
					}
				}
				sig2 = arma::exp(uobs * eta);

				/*********
				Sample rho
				*********/
				double xi = 0.5 * std::log((1.0 + rho) / (1.0 - rho));
				auto fx_rho = [&](double xi_input[])->double {
					return -loglik_rho(xi_input[0], wobs, xb, sig2);
				};
				double start[] = { xi };
				double xmin[] = { 0.0 };
				double ynewlo = 0.0;
				double reqmin = 1.0e-20;
				int konvge = 5;
				int kcount = 1000;
				double step[] = { 0.02 };
				int icount = 0;
				int numres = 0;
				int ifault = 0;
				nelmin(fx_rho, 1, start, xmin, &ynewlo, reqmin, step, konvge, kcount, &icount, &numres, &ifault);
				double xmax = xmin[0];
				double minll = ynewlo;

				mat cl(5,3, fill::zeros);
				vec dl(5, fill::zeros);
				double step_size = 0.2;
				rho_sampling_block:
					for (int iii=0; iii < 5; ++iii) {
						double e1 = static_cast<double>(iii-2);
						cl(iii,0) = std::pow(xmax + e1 * step_size, 2.0);
						cl(iii,1) = xmax + e1 * step_size;
						cl(iii,2) = 1.0;
						dl(iii) = -loglik_rho(xmax + e1 * step_size, wobs, xb, sig2);
					}

				for (int ni=0; ni < 5; ++ni) {
					if ((ni+1) != 3) {
						if (dl(ni) <= minll) {
							step_size *= 1.2;
							goto rho_sampling_block;
						}
					}
				}

				vec fl = arma::solve(cl.t() * cl, cl.t() * dl);
				double sigmaa = std::sqrt(0.5 / fl(0));


				double xi_prop = ::norm_rand() * sigmaa + xmax;

				// log-likelihood difference
				double ll_diff = loglik_rho(xi_prop, wobs, xb, sig2) -
								 loglik_rho(xi, wobs, xb, sig2) -
								 0.5 * (std::pow(xi - xmax, 2.0) - std::pow(xi_prop - xmax, 2.0)) / std::pow(sigmaa, 2.0);
				if (std::log(::unif_rand()) < ll_diff) {
					rho = (std::exp(2.0 * xi_prop) - 1.0) / (std::exp(2.0 * xi_prop) + 1.0);
					++rho_rates;
				}
			}
			beta_save.col(ikeep) = beta;
			alpha_save.col(ikeep) = alpha;
			eta_save.col(ikeep) = eta;
			lam_save.col(ikeep) = lam;
			sig2_save.col(ikeep) = sig2;
			rho_save(ikeep) = rho;
			prog.increment();
		}
	}
	alpha_rates /= static_cast<double>(ndiscard+nkeep*nskip);
	eta_rates /= static_cast<double>(ndiscard+nkeep*nskip);
	rho_rates /= static_cast<double>(ndiscard+nkeep*nskip);
	
	return ListBuilder()
	.add("beta", beta_save)
	.add("alpha", alpha_save)
	.add("eta", eta_save)
	.add("sig2", sig2_save)
	.add("rho", rho_save)
	.add("lam", lam_save)
	.add("rho_acceptance", rho_rates)
	.add("alpha_acceptance", alpha_rates)
	.add("eta_acceptance", eta_rates);
}