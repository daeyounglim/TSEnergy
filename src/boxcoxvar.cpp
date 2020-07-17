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
// [[Rcpp::depends(RcppArmadillo,RcppProgress))]]

Rcpp::List boxcoxvar(const arma::mat &yobs,
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
    const int Kp = K * p;

    // initialize parameters
    mat ystar = yobs;
    vec beta(Kp, fill::zeros);
    mat B = arma::reshape(beta, K, p);
    mat I_K(K, K, fill::eye);
    cube xtilde(K, Kp, T, fill::zeros);
    for (int t = 0; t < T; ++t) {
        rowvec xt = xobs.row(t);
        xtilde.slice(t) = arma::kron(xt, I_K);
    }
    mat A(K, K);
    std::generate(A.begin(), A.end(), ::unif_rand);
    mat C(K, p); C.fill(0.01);
    mat lam = C * xobs.t(); // K * T matrix
    
    mat Sig(K, K, fill::eye);
    mat Siginv = Sig / 5.0;
    Sig *= 5.0;
    
    mat w(arma::size(yobs), fill::zeros); // K * T matrix
    // bool miss_flag = arma::any(arma::any(miss));
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

    // Miscellaneous
    mat Bx = B * xobs.t();
    const double sigb2_i = 1.0 / sigb2;

    /*************************
	Parameters for adaptive MH
	*************************/
    const double obj_rate = 0.44;
    const int batch_length = std::min(50, ndiscard + nkeep * nskip);
    const int batch_total = (ndiscard + nskip * nkeep) / batch_length;

    mat A_accepts(K, K, fill::zeros);
    cube A_rates(K, K, batch_total, fill::zeros);
    mat A_tunings(K, K); A_tunings.fill(0.1);
    int batch_num = 0;

    mat C_accepts(K, p, fill::zeros);
    cube C_rates(K, p, batch_total, fill::zeros);
    mat C_tunings(K, p); C_tunings.fill(0.1);

    /************
	Miscellaneous
	************/
    mat xtilde1 = xtilde.slice(0);

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
            Rcpp::Rcout << idiscard + 1 << std::endl;
            /**************
             Update beta
             (beta =vec(B))
            **************/
            mat Sigbinv_post = xtilde1.t() * Siginv * xtilde1 +  eye<mat>(Kp, Kp) * sigb2_i;
            vec mub_post = xtilde1.t() * Siginv * w.col(0);
            for (int t = 1; t < T; ++t) {
                mat xtilde_t = xtilde.slice(t);
                Sigbinv_post += xtilde_t.t() * Siginv * xtilde_t;
                mub_post += xtilde_t.t() * Siginv * (w.col(t) - A * w.col(t-1));
            }
            mat Sigbinv_chol = arma::chol(Sigbinv_post);
            vec mub = arma::solve(arma::trimatu(Sigbinv_chol), arma::solve(trimatl(Sigbinv_chol.t()), mub_post));
            vec btmp(Kp);
            std::generate(btmp.begin(), btmp.end(), ::norm_rand);
            beta = mub + arma::solve(arma::trimatu(Sigbinv_chol), btmp);
            B = arma::reshape(beta, K, p);
            Bx = B * xobs.t();
            mat resid = w - Bx;

            /*******
            Update A
            ********/
            // double loglik_orig = loglik_A(A, w, Bx, Siginv);
            for (int k1 = 0; k1 < K; ++k1) {
                for (int k2 = 0; k2 < K; ++k2) {
                    double zij = 0.5 * std::log((1.0 + A(k1,k2)) / (1.0 - A(k1,k2)));
                    double zij_prop = ::norm_rand() * std::exp(A_tunings(k1,k2)) + zij;

                    // log-likelihood difference
                    double ll_diff = loglik_A_prop(zij_prop, k1, k2, A, w, Bx, Siginv) - loglik_A_prop(zij, k1, k2, A, w, Bx, Siginv);
                    if (std::log(::unif_rand()) < ll_diff) {
                        A(k1,k2) = (std::exp(2.0 * zij_prop) - 1.0) / (std::exp(2.0 * zij_prop) + 1.0);
                        // loglik_orig = loglik_A(A, w, Bx, Siginv); // recalculate the original log-likelihood
                        ++A_accepts(k1,k2);
                    }
                }
            }

            /*******
            Update C
            ********/
            double log_jacobian = arma::dot(lam - 1.0, yobs - a);
            double loglik_orig = loglik_C(lam, resid, w, A, Siginv, log_jacobian);
            for (int k = 0; k < K; ++k) {
                for (int j = 0; j < p; ++j) {
                    double cij = C(k,j);
                    double cij_prop = ::norm_rand() * std::exp(C_tunings(k,j)) + cij;

                    // log-likelihood difference
                    double ll_diff = loglik_C_prop(cij_prop, k, j, C, A, yobs, ystar, xobs, Bx, Siginv, a, sigc2) - loglik_orig - R::dnorm4(cij, 0.0, std::sqrt(sigc2), 1);
                    if (std::log(::unif_rand()) < ll_diff) {
                        C(k,j) = cij_prop;
                        lam.row(k) = C.row(k) * xobs.t(); // only kth row affected by the update
                        // lam = C * xobs.t();
                        for (int t = 0; t < T; ++t) {
                            if (lam(k,t) == 0.0) {
                                w(k,t) = std::log(ystar(k,t) - a);
                            } else {
                                w(k,t) = (std::pow(ystar(k,t) - a, lam(k,t)) - 1.0) / lam(k,t);
                            }
                        }
                        resid = w - Bx;
                        log_jacobian = arma::dot(lam - 1.0, yobs - a); // recalculate log of Jacobian
                        loglik_orig = loglik_C(lam, resid, w, A, Siginv, log_jacobian); // recalculate the original log-likelihood
                        ++C_accepts(k,j);
                    }
                }
            }

            /***********
            Update Sigma
            ***********/
            vec xx = resid.col(0);
            mat S = xx * xx.t();
            for (int t = 1; t < T; ++t) {
                xx = resid.col(t) - A * w.col(t-1);
                S += xx * xx.t();
            }
            Sig = riwish(static_cast<double>(T) + 1.0, S);
            Siginv = Sig.i();

            if ((idiscard+1) % batch_length == 0) {
                A_accepts /= static_cast<double>(batch_length);
                A_rates.slice(batch_num) = A_accepts;

                for (int k1 = 0; k1 < K; ++k1) {
                    for (int k2 = 0; k2 < K; ++k2) {
                        if (A_tunings(k1, k2) > obj_rate) {
                            A_tunings(k1, k2) += std::min(0.01, 1.0 / std::sqrt(static_cast<double>(batch_num+1)));
                        } else {
                            A_tunings(k1, k2) -= std::min(0.01, 1.0 / std::sqrt(static_cast<double>(batch_num+1)));
                        }
                    }
                }

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
                A_accepts.fill(0.0);
                C_accepts.fill(0.0);
            }
            prog.increment();
        }
    }

    cube A_save(K, K, nkeep, fill::zeros);
    cube B_save(K, p, nkeep, fill::zeros);
    cube C_save(K, p, nkeep, fill::zeros);
    cube lam_save(K, T, nkeep, fill::zeros);
    cube Sig_save(K, K, nkeep, fill::zeros);

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
                /**************
                 Update beta
                 (beta =vec(B))
                **************/
                mat Sigbinv_post = xtilde1.t() * Siginv * xtilde1 +  eye<mat>(Kp, Kp) * sigb2_i;
                vec mub_post = xtilde1.t() * Siginv * w.col(0);
                for (int t = 1; t < T; ++t) {
                    mat xtilde_t = xtilde.slice(t);
                    Sigbinv_post += xtilde_t.t() * Siginv * xtilde_t;
                    mub_post += xtilde_t.t() * Siginv * (w.col(t) - A * w.col(t-1));
                }
                mat Sigbinv_chol = arma::chol(Sigbinv_post);
                vec mub = arma::solve(arma::trimatu(Sigbinv_chol), arma::solve(trimatl(Sigbinv_chol.t()), mub_post));
                vec btmp(Kp);
                std::generate(btmp.begin(), btmp.end(), ::norm_rand);
                beta = mub + arma::solve(arma::trimatu(Sigbinv_chol), btmp);
                B = arma::reshape(beta, K, p);
                Bx = B * xobs.t();
                mat resid = w - Bx;

                /*******
                Update A
                ********/
                // double loglik_orig = loglik_A(A, w, Bx, Siginv);
                for (int k1 = 0; k1 < K; ++k1) {
                    for (int k2 = 0; k2 < K; ++k2) {
                        double zij = 0.5 * std::log((1.0 + A(k1,k2)) / (1.0 - A(k1,k2)));
                        double zij_prop = ::norm_rand() * std::exp(A_tunings(k1,k2)) + zij;

                        // log-likelihood difference
                        double ll_diff = loglik_A_prop(zij_prop, k1, k2, A, w, Bx, Siginv) - loglik_A_prop(zij, k1, k2, A, w, Bx, Siginv);
                        if (std::log(::unif_rand()) < ll_diff) {
                            A(k1,k2) = (std::exp(2.0 * zij_prop) - 1.0) / (std::exp(2.0 * zij_prop) + 1.0);
                            // loglik_orig = loglik_A(A, w, Bx, Siginv); // recalculate the original log-likelihood
                            ++A_accepts(k1,k2);
                        }
                    }
                }

                /*******
                Update C
                ********/
                double log_jacobian = arma::dot(lam - 1.0, yobs - a);
                double loglik_orig = loglik_C(lam, resid, w, A, Siginv, log_jacobian);
                for (int k = 0; k < K; ++k) {
                    for (int j = 0; j < p; ++j) {
                        double cij = C(k,j);
                        double cij_prop = ::norm_rand() * std::exp(C_tunings(k,j)) + cij;

                        // log-likelihood difference
                        double ll_diff = loglik_C_prop(cij_prop, k, j, C, A, yobs, ystar, xobs, Bx, Siginv, a, sigb2) - loglik_orig - R::dnorm4(cij, 0.0, std::sqrt(sigb2), 1);
                        if (std::log(::unif_rand()) < ll_diff) {
                            C(k,j) = cij_prop;
                            lam.row(k) = C.row(k) * xobs.t(); // only kth row affected by the update
                            // lam = C * xobs.t();
                            for (int t = 0; t < T; ++t) {
                                if (lam(k,t) == 0.0) {
                                    w(k,t) = std::log(ystar(k,t) - a);
                                } else {
                                    w(k,t) = (std::pow(ystar(k,t) - a, lam(k,t)) - 1.0) / lam(k,t);
                                }
                            }
                            resid = w - Bx;
                            log_jacobian = arma::dot(lam - 1.0, yobs - a); // recalculate log of Jacobian
                            loglik_orig = loglik_C(lam, resid, w, A, Siginv, log_jacobian); // recalculate the original log-likelihood
                            ++C_accepts(k,j);
                        }
                    }
                }

                /***********
                Update Sigma
                ***********/
                vec xx = resid.col(0);
                mat S = xx * xx.t();
                for (int t = 1; t < T; ++t) {
                    xx = resid.col(t) - A * w.col(t-1);
                    S += xx * xx.t();
                }
                Sig = riwish(static_cast<double>(T) + 1.0, S);
                Siginv = Sig.i();

                if (icount_mh % batch_length == 0) {
                    A_accepts /= static_cast<double>(batch_length);
                    A_rates.slice(batch_num) = A_accepts;

                    for (int k1 = 0; k1 < K; ++k1) {
                        for (int k2 = 0; k2 < K; ++k2) {
                            if (A_tunings(k1, k2) > obj_rate) {
                                A_tunings(k1, k2) += std::min(0.01, 1.0 / std::sqrt(static_cast<double>(batch_num+1)));
                            } else {
                                A_tunings(k1, k2) -= std::min(0.01, 1.0 / std::sqrt(static_cast<double>(batch_num+1)));
                            }
                        }
                    }

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
                    A_accepts.fill(0.0);
                    C_accepts.fill(0.0);
                }
            }

            A_save.slice(ikeep) = A;
            B_save.slice(ikeep) = B;
            C_save.slice(ikeep) = C;
            lam_save.slice(ikeep) = lam;
            Sig_save.slice(ikeep) = Sig;
            prog.increment();
        }
    }

    return ListBuilder()
    .add("A", A_save)
    .add("B", B_save)
    .add("C", C_save)
    .add("Sig", Sig_save)
    .add("lam", lam_save)
    .add("A_acceptance", A_rates)
    .add("C_acceptance", C_rates);
}