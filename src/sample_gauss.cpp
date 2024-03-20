#define RCPPDIST_DONT_USE_ARMA

#include <RcppEigen.h>
#include <Rcpp.h>
#include <complex>
#include <RcppDist.h>
#include "prior.h"
#include "MCMCtools.h"

extern "C"
{
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_hyperg.h>
}
using namespace Rcpp;
using namespace Eigen;

double l_lik(Eigen::VectorXd res,
             double s2e){
  double out =  - (res.size() / 2.0) * std::log(s2e) - ((1.0 / (2.0 * s2e)) * res.squaredNorm());
  return(out);
}


void gibbs_beta(Eigen::VectorXd& beta,
                Eigen::VectorXd& Xbeta,
                Eigen::MatrixXd XtX,
                Eigen::MatrixXd X,
                Eigen::MatrixXd Q_beta,
                double s2e,
                Eigen::VectorXd res_beta){

  Eigen::LLT<MatrixXd> chol_Q_fc(XtX / s2e +  Q_beta);
  Eigen::VectorXd mu = chol_Q_fc.solve((1.0 / s2e) * X.transpose() * res_beta);
  Eigen::VectorXd z = Rcpp::as<VectorXd>(Rcpp::rnorm(X.cols()));

  beta = chol_Q_fc.matrixU().solve(z) + mu;
  Xbeta = X * beta;

}


void metropolis_s2e(double& s2e,
                    Eigen::VectorXd res_eps,
                    Rcpp::List pri_s2e,
                    double FFe,
                    int& n_acc_s2e){
  double s2e_prop;

  double l_fc_s2e_curr = l_lik(res_eps, s2e) + l_pri_s2(s2e, pri_s2e);
  rprop_rw(s2e_prop, s2e, FFe);
  double l_fc_s2e_prop = l_lik(res_eps, s2e_prop) + l_pri_s2(s2e_prop, pri_s2e);

  if(metropolis_acceptance(l_fc_s2e_prop, l_fc_s2e_curr) == 1L){
    s2e = s2e_prop;
    n_acc_s2e +=1;
  }

}


void compute_Qb_g_fc(Eigen::SparseMatrix<double>& Q,
                     Eigen::VectorXd& b,
                     Eigen::SparseMatrix<double> ZtZ,
                     Eigen::SparseMatrix<double> K_g,
                     Eigen::SparseMatrix<double> Z,
                     Eigen::VectorXd res_g,
                     double s2e, double s2g){

  Q = ZtZ/s2e + K_g/s2g;
  b = (1.0 / s2e) * Z.transpose() * res_g;

}


double lfc_block(Eigen::VectorXd gamma,
                 double s2g,
                 Eigen::VectorXd mu,
                 Eigen::SparseMatrix<double> Q,
                 double l_det_Q,
                 Eigen::MatrixXd A,
                 Eigen::VectorXd e,
                 Eigen::MatrixXd W,
                 Eigen::SparseMatrix<double> K,
                 int rank,
                 Rcpp::List pri_s2,
                 double s2e,
                 Eigen::VectorXd res_g,
                 Eigen::SparseMatrix<double> Z){

  double lfc_gamma, lfc_joint, lfc;

  Eigen::VectorXd g_mu = gamma - mu;

  lfc_gamma = 0.5 * l_det_Q - 0.5 * g_mu.transpose() * Q * g_mu;

  if (!NumericVector::is_na(A(0,0))){
    Eigen::VectorXd e_Amu = e - A * mu;
    lfc_gamma -= - 0.5 * std::log(W.determinant()) - 0.5 * e_Amu.transpose() * W.inverse() * e_Amu;
  }

  lfc_joint = l_lik(res_g - Z * gamma, s2e) - (rank / 2.0) * std::log(s2g) -
    (1.0 / (2.0 * s2g)) * gamma.transpose() * K * gamma +
    l_pri_s2(s2g, pri_s2);
  lfc = lfc_joint - lfc_gamma;

  return(lfc);
}

//'@export
// [[Rcpp::export]]
Rcpp::List sample_blmm(const Eigen::VectorXd y,
                       const Eigen::MatrixXd X,
                       const Rcpp::List Z_list,
                       const Rcpp::List K_list,
                       const std::vector<int> rank_K_g,
                       const Rcpp::List A_list,
                       const Rcpp::List e_list,
                       const Eigen::MatrixXd S_beta,
                       const Rcpp::List pri_s2e,
                       const Rcpp::List pri_s2b,
                       const Rcpp::List pri_s2g,
                       const Eigen::VectorXd beta_init,
                       const Rcpp::List g_init_list,
                       const Eigen::VectorXd S2g_init,
                       const double s2e_init,
                       double FFe,
                       Eigen::VectorXd FFb,
                       Eigen::VectorXd FFg,
                       const int niter, const int pr, const int thin,
                       const int ntuning, const int stop_tuning){

  ///////////////////////////////////////////////////////////////////////
  // Initialise mcmc objects - s2e, s2b and s2g
  ///////////////////////////////////////////////////////////////////////

  ////////////// s2e  //////////////
  double s2e = s2e_init;

  int
    Ta_s2e = 0, n_acc_s2e = 0;

  ////////////// s2b  //////////////
  int p = X.cols();

  Eigen::VectorXd
    s2b = S_beta.diagonal();

  std::vector<int>
    Ta_s2b(p, 0), n_acc_s2b(p, 0);

  ////////////// s2g  //////////////
  int q = Z_list.length();

  double
    s2g_prop;

  Eigen::VectorXd
    s2g_curr = S2g_init;

  std::vector<int>
    n_acc_s2g(q);

  std::vector<int>
    Ta_s2g(q);


  ///////////////////////////////////////////////////////////////////////
  // Initialise mcmc objects - beta
  ///////////////////////////////////////////////////////////////////////
  Eigen::VectorXd
    beta = beta_init, Xbeta = X * beta;

  Eigen::MatrixXd
    Q_beta = S_beta.inverse();

  //////////////////////////////////////////////////////////////////////
  // Fixed quantities - beta
  ///////////////////////////////////////////////////////////////////////
  Eigen::MatrixXd
    XtX = X.transpose() * X;

  ///////////////////////////////////////////////////////////////////////
  // Initialise mcmc objects - gamma
  ///////////////////////////////////////////////////////////////////////
  std::vector<Eigen::VectorXd>
    g_curr(q);

  for(int j=0; j<q ;j++){
    g_curr[j] = Rcpp::as<Eigen::VectorXd>(g_init_list[j]);
  }

  int n = X.rows();

  Eigen::VectorXd
    Zg = Eigen::VectorXd::Zero(n),
      g_prop, res_g, b_g,
      mu_g_fc;

  Eigen::SparseMatrix<double>
    Q_g_fc;

  Eigen::MatrixXd
    V, W;

  double
    ldet_Q_g_fc,
    l_fc_g_s2g_prop, l_fc_g_s2g_curr;

  ///////////////////////////////////////////////////////////////////////
  // Initialise fixed quantities - gamma
  ///////////////////////////////////////////////////////////////////////
  std::vector<Eigen::SparseMatrix<double> >
    Z(q), ZtZ(q), K_g(q);

  std::vector<Eigen::MatrixXd>
    A(q);

  std::vector<Eigen::VectorXd>
    e(q);

  std::vector<int>
    mj(q);

  for(int j=0; j<q ;j++){
    Z[j] = Rcpp::as<Eigen::SparseMatrix<double> >(Z_list[j]);
    ZtZ[j] = Z[j].transpose()*Z[j];
    mj[j] = Z[j].cols();
    K_g[j] = Rcpp::as<Eigen::SparseMatrix<double> >(K_list[j]);
    A[j] = Rcpp::as<Eigen::MatrixXd>(A_list[j]);
    e[j] = Rcpp::as<Eigen::VectorXd>(e_list[j]);
  }

  std::vector<Eigen::SimplicialLLT<Eigen::SparseMatrix<double>, Eigen::Lower> >
    Chol_Q_g(q);

  for(int j=0; j<q ;j++){
    compute_res_gamma(res_g, j, y, Xbeta, Z, g_curr);
    compute_Qb_g_fc(Q_g_fc, b_g, ZtZ[j], K_g[j], Z[j], res_g, s2e, s2g_curr[j]);
    Chol_Q_g[j].analyzePattern(Q_g_fc);
  }


  ///////////////////////////////////////////////////////////////////////
  // Initialise objects for storage
  ///////////////////////////////////////////////////////////////////////
  int storeind = 0, nstore = niter / thin;

  Eigen::VectorXd
    s2e_mc(nstore);

  Eigen::MatrixXd
    beta_mc(nstore, p),
    s2b_mc(nstore, p);

  std::vector<double>
    FFes;

  Eigen::MatrixXd
    FFbs(stop_tuning / ntuning, p);

  Eigen::MatrixXd
    s2g_mc(nstore, q);

  std::vector<Eigen::MatrixXd>
    g_mc(q);
  for(int j=0; j<q; j++){
    g_mc[j] = Eigen::MatrixXd::Zero(nstore, mj[j]);
  }

  Eigen::MatrixXd
    FFgs(nstore, q);


  ///////////////////////////////////////////////////////////////////////
  // MCMC sampler
  ///////////////////////////////////////////////////////////////////////
  time_t now;
  for(int k=0;k<niter;k++){

    ///////////////////////////////////////////////////////////////////////
    // Sampling beta and s2b
    ///////////////////////////////////////////////////////////////////////
    gibbs_beta(beta, Xbeta, XtX, X, Q_beta, s2e, y - Zg);
    metropolis_s2b(s2b, Q_beta, beta, pri_s2b, FFb, n_acc_s2b);

    ///////////////////////////////////////////////////////////////////////
    // Sampling (gamma, s2g)
    ///////////////////////////////////////////////////////////////////////
    Zg = Zg * 0.0;
    for(int j=0; j<q ;j++){
      compute_res_gamma(res_g, j, y, Xbeta, Z, g_curr);

      rprop_rw(s2g_prop, s2g_curr[j], FFg[j]);

      // Full conditionals parameters - curr
      compute_Qb_g_fc(Q_g_fc, b_g, ZtZ[j], K_g[j], Z[j], res_g, s2e, s2g_curr[j]);
      Chol_Q_g[j].factorize(Q_g_fc);
      mu_g_fc = Chol_Q_g[j].solve(b_g);
      ldet_Q_g_fc = 2.0 * Chol_Q_g[j].matrixL().toDense().diagonal().array().log().sum();
      if (!Rcpp::NumericVector::is_na(A[j](0,0))){
        V = Chol_Q_g[j].solve(A[j].transpose());
        W = A[j] * V;
        //g_curr[j] è già costretto
      }
      l_fc_g_s2g_curr = lfc_block(g_curr[j], s2g_curr[j],
                                  mu_g_fc, Q_g_fc, ldet_Q_g_fc, A[j], e[j], W,
                                  K_g[j], rank_K_g[j], pri_s2g[j], s2e, res_g, Z[j]);

      // Full conditionals parameters - prop
      compute_Qb_g_fc(Q_g_fc, b_g, ZtZ[j], K_g[j], Z[j], res_g, s2e, s2g_prop);
      Chol_Q_g[j].factorize(Q_g_fc);
      mu_g_fc = Chol_Q_g[j].solve(b_g);
      ldet_Q_g_fc = 2.0 * Chol_Q_g[j].matrixL().toDense().diagonal().array().log().sum();

      // proposal of gamma
      g_prop = Chol_Q_g[j].permutationPinv() *
        Chol_Q_g[j].matrixU().solve(Rcpp::as<Eigen::VectorXd>(Rcpp::rnorm(K_g[j].cols())));
      g_prop += mu_g_fc;

      if (!Rcpp::NumericVector::is_na(A[j](0,0))){
        V = Chol_Q_g[j].solve(A[j].transpose());
        W = A[j] * V;
        correction_sample_ulc(g_prop, A[j], e[j], W, V);
      }
      l_fc_g_s2g_prop = lfc_block(g_prop, s2g_prop,
                                  mu_g_fc, Q_g_fc, ldet_Q_g_fc, A[j], e[j], W,
                                  K_g[j], rank_K_g[j], pri_s2g[j], s2e, res_g, Z[j]);

      // acceptance step
      if(metropolis_acceptance(l_fc_g_s2g_prop, l_fc_g_s2g_curr) == 1L){
        s2g_curr[j] = s2g_prop;
        g_curr[j] = g_prop;
        n_acc_s2g[j] += 1;
      }

      Zg += Z[j]*g_curr[j];
    }

    ///////////////////////////////////////////////////////////////////////
    // Sampling s2e
    ///////////////////////////////////////////////////////////////////////
    metropolis_s2e(s2e, y - Xbeta - Zg, pri_s2e, FFe, n_acc_s2e);

    ///////////////////////////////////////////////////////////////////////
    // Tuning
    ///////////////////////////////////////////////////////////////////////
    if(k < stop_tuning){
      if((k+1)%ntuning==0) {
        tune(FFe, n_acc_s2e, ntuning, 0.44, Ta_s2e);
        FFes.push_back(FFe);
        for(int j = 0; j<p ;j++){
          if(is_mixture_pri_beta(pri_s2b[j]) == 1L){
            tune(FFb[j], n_acc_s2b[j], ntuning, 0.44, Ta_s2b[j]);
            FFbs(Ta_s2b[j]-1, j) = FFb[j];
          }
        }
        for(int j=0; j<q ;j++){
          tune(FFg[j], n_acc_s2g[j], ntuning, 0.44, Ta_s2g[j]);
          FFgs(Ta_s2g[j]-1, j) = FFg[j];
        }
      }
    }

    ///////////////////////////////////////////////////////////////////////
    // Storage
    ///////////////////////////////////////////////////////////////////////
    if ((k+1)%thin == 0){
      s2e_mc(storeind) = s2e;
      beta_mc.row(storeind) = beta;
      s2b_mc.row(storeind) = s2b;
      s2g_mc.row(storeind) = s2g_curr;
      for(int j=0;j<q;j++){
        g_mc[j].row(storeind) = g_curr[j];
      }
      storeind += 1;
    }

    if ((k+1)%pr==0){
      time(&now);
      Rprintf("Iteration: %d: %s", k+1, ctime(&now));
    }

  }

  return Rcpp::List::create(
    Rcpp::Named("s2e") = s2e_mc,
    Rcpp::Named("s2b") = s2b_mc,
    Rcpp::Named("s2g") = s2g_mc,
    Rcpp::Named("beta") = beta_mc,
    Rcpp::Named("gamma") = g_mc,
    Rcpp::Named("FFes") = FFes,
    Rcpp::Named("FFbs") = FFbs,
    Rcpp::Named("FFgs") = FFgs,
    Rcpp::Named("Ta_s2b") = Ta_s2b,
    Rcpp::Named("Qbeta") = Q_beta);
}

