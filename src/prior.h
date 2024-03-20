#ifndef PRIOR_H
#define PRIOR_H

#include <RcppEigen.h>
#include <Rcpp.h>
#include <complex>
#include <RcppDist.h>
#include "Hypergeo_fun.h"

extern "C"
{
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_hyperg.h>
}
using namespace Rcpp;
using namespace Eigen;

double l_pri_s2(const double s, Rcpp::List pri_s2);


#endif
