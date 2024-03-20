#ifndef HYPERGEO_FUN_H
#define HYPERGEO_FUN_H

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

Eigen::VectorXd l_2F1(double q_pri, double p_pri,
                      double alpha, double alpha_til,
                    double s, double x, double lmlt);


#endif
