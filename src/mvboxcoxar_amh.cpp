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
#include "misc_mvboxcoxar.h"
#include "nelmin.h"
// [[Rcpp::depends(RcppArmadillo,RcppProgress))]]

Rcpp::List mvboxcoxar_amh(const arma::mat &yobs,
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
    mat A(K, K, fill::zeros);
    mat C(K, p); C.fill(0.01);
    mat lam = C * xobs.t(); // K * T matrix
    
    mat Sig(K, K, fill::eye);
    mat Siginv = Sig;
    
    mat w(arma::size(yobs), fill::zeros); // K * T matrix
    bool miss_flag = any(miss);
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
    const double sigc2_i = 1.0 / sigc2;

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

            /**************
             Update beta
             (beta =vec(B))
            **************/
            mat Sigbinv_post = xtilde1 * Siginv * xtilde1.t() +  eye<mat>(Kp, Kp) * sigb2_i;
            vec mub_post = xtilde1 * Siginv * w.col(0);
            for (int t = 1; t < T; ++t) {
                mat xtilde_t = xtilde.slice(t);
                Sigbinv_post += xtilde_t * Siginv * xtilde_t.t();
                mub_post += xtilde_t * Siginv * (w.col(t) - A * w.col(t-1));
            }
            mat Sigbinv_chol = arma::chol(Sigbinv_post);
            vec mub = arma::solve(arma::trimatu(Sigbinv_chol), arma::solve(trimatl(Sigbinv_chol.t()), mub_post));
            vec btmp(Kp);
            std::generate(btmp.begin(), btmp.end(), ::norm_rand);
            beta = mub + arma::solve(arma::trimatu(Sigbinv_chol), btmp);
            B = arma::reshape(beta, K, p);
            Bx = B * xobs.t();

            /*******
            Update A
            ********/
            for (int k1 = 0; k1 < K; ++k1) {
                for (int k2 = 0; k2 < K; ++k2) {
                    
                }
            }
        }
    }