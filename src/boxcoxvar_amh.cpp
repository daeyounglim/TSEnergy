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
#include "misc_mvboxcoxar.h"
#include "nelmin.h"
// [[Rcpp::depends(RcppArmadillo,RcppProgress))]]


Rcpp::List boxcoxvar_amh(const arma::mat &yobs,
                         const arma::umat &miss,
                         const double &a,
                         const arma::mat &xobs,
                         const double &sigb2,
                         const double &sigc2,
                         const int &ndiscard,
                         const int &nskip,
                         const int &nkeep,
                         const bool verbose)
{
    using namespace arma;
    using namespace Rcpp;
    using namespace R;
    using namespace std;

    const int T = yobs.n_cols; // yobs is a K * T matrix
    const int K = yobs.n_rows;
    const int p = xobs.n_cols; // xobs is a T * p matrix
    const int sizeb = K * (K+p);

    // initialize parameters
    mat ystar = yobs;
    vec b(sizeb, fill::zeros);
    mat I_K(K, K, fill::eye);
    mat xtilde(T, p+K, fill::zeros);
    mat C(K, p); C.fill(0.01);
    mat lam = C * xobs.t(); // K * T matrix

    mat Sig(K, K, fill::eye);
    mat Siginv = Sig / 5.0;
    Sig *= 5.0;

    mat w(arma::size(yobs), fill::zeros); // K * T matrix
    bool miss_flag = arma::any(arma::any(miss));
    for (int t = 0; t < T; ++t) {
        for (int k = 0; k < K; ++k) {
            if (miss(k,t)) {
                double normrv = ::norm_rand();
                ystar(k,t) = normrv;
            }
            if (lam(k,t) == 0.0) {
                w(k,t) = std::log(ystar(k,t) - a);
            } else {
                w(k,t) = (std::pow(ystar(k,t) - a, lam(k,t)) - 1.0) / lam(k,t);
            }
        }
    }
    cube Xts(K, sizeb, T, fill::zeros);
    for (int t = 0; t < T; ++t) {
    	rowvec xtilde_i(p+K, fill::zeros);
    	xtilde_i.head(p) = xobs.row(t);
    	xtilde_i.tail(K) = arma::trans(w.col(t));

    	xtilde.row(t) = xtilde_i;
    	Xts.slice(t) = arma::kron(xtilde_i, I_K);
    }

    /*************************
	Parameters for adaptive MH
	*************************/
    const double obj_rate = 0.44;
    const int batch_length = std::min(50, ndiscard + nkeep * nskip);
    const int batch_total = (ndiscard + nskip * nkeep) / batch_length;

    mat C_accepts(K, p, fill::zeros);
    cube C_rates(K, p, batch_total, fill::zeros);
    mat C_tunings(K, p); C_tunings.fill(0.1);

    /*******************
	Begin burn-in period
	*******************/
    if (verbose)
    {
        Rcout << "Warming up" << endl;
    }
    {
        Progress prog(ndiscard, verbose);
        for (int idiscard = 0; idiscard < ndiscard; ++idiscard)
        {
            if (Progress::check_abort())
            {
                return Rcpp::List::create(Rcpp::Named("error") = "user interrupt aborted");
            }
            /*******
            Update b
            *******/
            mat Sigb = eye<mat>(sizeb, sizeb) / sigb2;
            vec mub_post(sizeb, fill::zeros);
            for (int t=0; t < T; ++t) {
            	mat Xt = Xts.slice(t);
            	Sigb += Xt.t() * Siginv * Xt;
            	mub_post += Xt.t() * Siginv * w.col(t);
            }
            mat Sigbinv_chol = arma::chol(Sigb);
            vec mub = arma::solve(arma::trimatu(Sigbinv_chol), arma::solve(trimatl(Sigbinv_chol.t()), mub_post));
            vec btmp(sizeb);
            std::generate(btmp.begin(), btmp.end(), ::norm_rand);
            b = mub + arma::solve(arma::trimatu(Sigbinv_chol), btmp);

            /*********
			Update Sig
            *********/
            mat S(K, K, fill::zeros);
            for (int t = 0; t < T; ++t) {
            	mat Xt = Xts.slice(t);
            	vec resid_t = w.col(t) - Xt * b;
            	S += resid_t * resid_t.t();
            }
            Siginv = riwish(static_cast<double>(T+1), S);

            /*******
            Update C
            *******/
            double log_jacobian = arma::dot(lam - 1.0, ystar - a);
            double loglik_orig = loglik_C(resid, Siginv, log_jacobian);
            for (int k = 0; k < K; ++k) {
                for (int j = 0; j < p; ++j) {
                    double cij = C(k,j);
                    double cij_prop = ::norm_rand() * std::exp(C_tunings(k,j)) + cij;

                    mat C_prop = C;
                    C_prop(k,j) = cij_prop;
                    mat lam_prop = lam;
                    lam_prop.row(k) = C_prop.row(k) * xobs.t(); // only kth row affected by the update

                    mat w_prop = w;
                    for (int t = 0; t < T; ++t) {
                        if (lam_prop(k,t) == 0.0) {
                            w_prop(k,t) = std::log(ystar(k,t) - a);
                        } else {
                            w_prop(k,t) = (std::pow(ystar(k,t) - a, lam_prop(k,t)) - 1.0) / lam_prop(k,t);
                        }
                    }

                    mat Xb_prop(K, T, fill::zeros);

                    for (int t = 0; t < T; ++t) {
				    	rowvec xtilde_i(p+K, fill::zeros);
				    	xtilde_i.head(p) = xobs.row(t);
				    	xtilde_i.tail(K) = arma::trans(w.col(t));

				    	xtilde.row(t) = xtilde_i;
				    	mat Xt = arma::kron(xtilde_i, I_K);
				    	Xb_prop.col(t) = Xt * b;
				    }


                    // log-likelihood difference
                    double ll_diff = loglik_C_prop_amh(C_prop, cij_prop,lam_prop, resid_prop, ystar, Siginv, a, sigc2) - 
                    				 loglik_orig - R::dnorm4(cij, 0.0, std::sqrt(sigc2), 1);
                    if (std::log(::unif_rand()) < ll_diff) {
                        C(k,j) = cij_prop;
                        lam = lam_prop;
                        w = w_prop;

                        resid = w - Xb_prop;
                        log_jacobian = arma::dot(lam - 1.0, ystar - a); // recalculate log of Jacobian
                        loglik_orig = loglik_C(resid, Siginv, log_jacobian); // recalculate the original log-likelihood
                        ++C_accepts(k,j);
                    }
                }
            }
            for (int t = 0; t < T; ++t) {
		    	rowvec xtilde_i(p+K, fill::zeros);
		    	xtilde_i.head(p) = xobs.row(t);
		    	xtilde_i.tail(K) = arma::trans(w.col(t));

		    	xtilde.row(t) = xtilde_i;
		    	Xts.slice(t) = arma::kron(xtilde_i, I_K);
		    }


            if ((idiscard+1) % batch_length == 0) {
                C_accepts /= static_cast<double>(batch_length);
                C_rates.slice(batch_num) = C_accepts;

                for (int i = 0; i < K; ++i) {
                    for (int j = 0; j < p; ++j) {
                        if (C_tunings(i,j) > obj_rate) {
                            C_tunings(i, j) += std::min(0.01, 1.0 / std::sqrt(static_cast<double>(batch_num+1)));
                        } else {
                            C_tunings(i, j) -= std::min(0.01, 1.0 / std::sqrt(static_cast<double>(batch_num+1)));
                        }
                    }
                }

                ++batch_num;
                C_accepts.fill(0.0);
            }
			prog.increment();
        }
    }


    /*******************
	Begin sampling period
	*******************/
	mat bsave(sizeb, nkeep, fill::zeros);
	cube Siginvsave(K, K, nkeep, fill::zeros);
	cube Csave(K, p, fill::zeros);

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

				/*******
	            Update b
	            *******/
	            mat Sigb = eye<mat>(sizeb, sizeb) / sigb2;
	            vec mub_post(sizeb, fill::zeros);
	            for (int t=0; t < T; ++t) {
	            	mat Xt = Xts.slice(t);
	            	Sigb += Xt.t() * Siginv * Xt;
	            	mub_post += Xt.t() * Siginv * w.col(t);
	            }
	            mat Sigbinv_chol = arma::chol(Sigb);
	            vec mub = arma::solve(arma::trimatu(Sigbinv_chol), arma::solve(trimatl(Sigbinv_chol.t()), mub_post));
	            vec btmp(sizeb);
	            std::generate(btmp.begin(), btmp.end(), ::norm_rand);
	            b = mub + arma::solve(arma::trimatu(Sigbinv_chol), btmp);

	            /*********
				Update Sig
	            *********/
	            mat S(K, K, fill::zeros);
	            for (int t = 0; t < T; ++t) {
	            	mat Xt = Xts.slice(t);
	            	vec resid_t = w.col(t) - Xt * b;
	            	S += resid_t * resid_t.t();
	            }
	            Siginv = riwish(static_cast<double>(T+1), S);

	            /*******
	            Update C
	            *******/
	            double log_jacobian = arma::dot(lam - 1.0, ystar - a);
	            double loglik_orig = loglik_C(resid, Siginv, log_jacobian);
	            for (int k = 0; k < K; ++k) {
	                for (int j = 0; j < p; ++j) {
	                    double cij = C(k,j);
	                    double cij_prop = ::norm_rand() * std::exp(C_tunings(k,j)) + cij;

	                    mat C_prop = C;
	                    C_prop(k,j) = cij_prop;
	                    mat lam_prop = lam;
	                    lam_prop.row(k) = C_prop.row(k) * xobs.t(); // only kth row affected by the update

	                    mat w_prop = w;
	                    for (int t = 0; t < T; ++t) {
	                        if (lam_prop(k,t) == 0.0) {
	                            w_prop(k,t) = std::log(ystar(k,t) - a);
	                        } else {
	                            w_prop(k,t) = (std::pow(ystar(k,t) - a, lam_prop(k,t)) - 1.0) / lam_prop(k,t);
	                        }
	                    }

	                    mat Xb_prop(K, T, fill::zeros);

	                    for (int t = 0; t < T; ++t) {
					    	rowvec xtilde_i(p+K, fill::zeros);
					    	xtilde_i.head(p) = xobs.row(t);
					    	xtilde_i.tail(K) = arma::trans(w.col(t));

					    	xtilde.row(t) = xtilde_i;
					    	mat Xt = arma::kron(xtilde_i, I_K);
					    	Xb_prop.col(t) = Xt * b;
					    }


	                    // log-likelihood difference
	                    double ll_diff = loglik_C_prop_amh(C_prop, cij_prop,lam_prop, resid_prop, ystar, Siginv, a, sigc2) - 
	                    				 loglik_orig - R::dnorm4(cij, 0.0, std::sqrt(sigc2), 1);
	                    if (std::log(::unif_rand()) < ll_diff) {
	                        C(k,j) = cij_prop;
	                        lam = lam_prop;
	                        w = w_prop;

	                        resid = w - Xb_prop;
	                        log_jacobian = arma::dot(lam - 1.0, ystar - a); // recalculate log of Jacobian
	                        loglik_orig = loglik_C(resid, Siginv, log_jacobian); // recalculate the original log-likelihood
	                        ++C_accepts(k,j);
	                    }
	                }
	            }
	            for (int t = 0; t < T; ++t) {
			    	rowvec xtilde_i(p+K, fill::zeros);
			    	xtilde_i.head(p) = xobs.row(t);
			    	xtilde_i.tail(K) = arma::trans(w.col(t));

			    	xtilde.row(t) = xtilde_i;
			    	Xts.slice(t) = arma::kron(xtilde_i, I_K);
			    }


	            if (icount_mh % batch_length == 0) {
                    C_accepts /= static_cast<double>(batch_length);
                    C_rates.slice(batch_num) = C_accepts;

                    for (int i = 0; i < K; ++i) {
                        for (int j = 0; j < p; ++j) {
                            if (C_tunings(i,j) > obj_rate) {
                                C_tunings(i, j) += std::min(0.01, 1.0 / std::sqrt(static_cast<double>(batch_num+1)));
                            } else {
                                C_tunings(i, j) -= std::min(0.01, 1.0 / std::sqrt(static_cast<double>(batch_num+1)));
                            }
                        }
                    }

                    ++batch_num;
                    C_accepts.fill(0.0);
                }
			}

			bsave.col(ikeep) = b;
			Siginvsave.slice(ikeep) = Siginv;
			Csave.slice(ikeep) = C;
		}
	}

	return ListBuilder()
	.add("b", bsave)
	.add("Siginv", Siginvsave)
	.add("C", Csave);
}
