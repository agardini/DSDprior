#ifndef MCMCTOOLS_H
#define MCMCTOOLS_H

#include <RcppEigen.h>
#include <Rcpp.h>
#include <complex>
#include <RcppDist.h>

extern "C"
{
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_hyperg.h>
}
using namespace Rcpp;
using namespace Eigen;

double my_rpg(double b, double c);

void rprop_rw(double& s2_prop,
              double s2_curr,
              double FF);


int metropolis_acceptance(double lfc_prop,
                          double lfc_curr);

void tune(double& FF,
          int& n_acc,
          int ntuning,
          double a_opt,
          int& Ta);


int is_mixture_pri_beta(Rcpp::List pri);

void metropolis_s2b(Eigen::VectorXd& s2b,
                    Eigen::MatrixXd& Q_beta,
                    Eigen::VectorXd beta,
                    Rcpp::List pri_s2b,
                    Eigen::VectorXd FFb,
                    std::vector<int>& n_acc_s2b);

void metropolis_s2g_noblock(Eigen::VectorXd& s2g,
                            Eigen::MatrixXd Q_gamma_j,
                            Eigen::VectorXd gamma_j,
                            Rcpp::List pri_s2g_j,
                            double rank_Q_j,
                            double FFg_j,
                            int j,
                            std::vector<int>& n_acc_s2g);

void compute_res_gamma(Eigen::VectorXd& res_gamma, int j,
                       Eigen::VectorXd y,
                       Eigen::VectorXd Xbeta,
                       std::vector<Eigen::SparseMatrix<double> > Z,
                       std::vector<Eigen::VectorXd> gamma_curr);

void correction_sample_ulc(Eigen::VectorXd& gamma,
                           Eigen::MatrixXd A,
                           Eigen::VectorXd e,
                           Eigen::MatrixXd W,
                           Eigen::MatrixXd V);





#endif
