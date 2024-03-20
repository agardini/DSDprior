#define RCPPDIST_DONT_USE_ARMA

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




bool is_neg_int(double x){
  if(x>0) {
    return(false);
  }
  if(std::floor(x) == x){
    return(true);
  }else{
    return(false);}
}


double lmolt_pri(double s, double alpha, double alpha_til, double beta, double beta_til,
                 double b_pri, double p_pri, double q_pri){

  double out = std::log(beta / (b_pri * beta_til)) +
    gsl_sf_lngamma(q_pri + p_pri) - gsl_sf_lngamma(p_pri) - gsl_sf_lngamma(q_pri) +
    gsl_sf_lngamma(q_pri + alpha) - gsl_sf_lngamma(alpha) +
    gsl_sf_lngamma(alpha_til) -  gsl_sf_lngamma(alpha_til + q_pri) +
    (-q_pri - 1.0) * std::log((s * beta) / (b_pri * beta_til));
  return(out);
}



double l_dbalpri_2f1(const double s,
                     const double alpha,
                     const double alpha_til,
                     const double beta,
                     const double beta_til,
                     const double b_pri,
                     const double p_pri,
                     const double q_pri){
  double x;
  double alpha_C = alpha;

  if(is_neg_int(p_pri - alpha)){
    alpha_C = alpha + alpha*0.0000000001;
  }
  if(is_neg_int(alpha_til - alpha)){
    alpha_C = alpha + alpha*0.0000000001;
  }


  x = beta / (b_pri * beta_til);

  double lmlt = lmolt_pri(s, alpha_C, alpha_til, beta, beta_til,
                          b_pri, p_pri, q_pri);
  Eigen::VectorXd l_f21 = l_2F1(q_pri, p_pri,
                                alpha_C, alpha_til,
                                s, x, lmlt);
  double out = lmlt + l_f21[1];
  //printf("s: %f; lval: %f", s, out);
  return(out);


}



//'@export
// [[Rcpp::export]]
double l_pri_s2(const double s, Rcpp::List pri_s2){
  double out;
  std::string check = pri_s2["Kind_pri"];

  if(check == "B2"){
    double b = pri_s2["b"], p = pri_s2["p"], q = pri_s2["q"];
    out = (p - 1.0) * std::log(s) - (p + q) * std::log(1.0 + s / b);
  };

  if(check == "InvGamma"){
    double a = pri_s2["a"], b = pri_s2["b"];
    out = -(a + 1.0) * std::log(s) - b / s;
  };

  if(check == "Gamma"){
    double a = pri_s2["a"], b = pri_s2["b"];
    out = (a - 1.0) * std::log(s) - b * s;
  };

  if(check == "Weibull"){
    double lambda = pri_s2["lambda"], kappa = pri_s2["kappa"];
    out = std::log(kappa / lambda) + (kappa - 1.0) *
      std::log(s / lambda) - std::pow(s / lambda, kappa);
  };

  if(check == "DSD"){
    double
    alpha = pri_s2["alpha"], alpha_til = pri_s2["alpha_til"],
                                               beta = pri_s2["beta"], beta_til = pri_s2["beta_til"],
                                                                                       b = pri_s2["b"], p = pri_s2["p"], q = pri_s2["q"];
    out = l_dbalpri_2f1(s, alpha, alpha_til, beta, beta_til,
                        b, p, q);
    //Rprintf("s: %f; lval: %f", s, out);
  };

  return(out);
}

