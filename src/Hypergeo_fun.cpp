#define RCPPDIST_DONT_USE_ARMA

#include <RcppEigen.h>
#include <RcppGSL.h>
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


int sign(double x) {
  return (x > 0) - (x < 0);
}

void f21_wrap(double a, double b, double c, double s, Eigen::VectorXd& res){
  gsl_set_error_handler_off();
  gsl_sf_result result;
  int status = 0;

  status = gsl_sf_hyperg_2F1_e(a, b, c, s, &result);

  res[0] = status;
  res[1] = result.val;
  res[2] = result.err;
  res[3] = 1.0;
}


double sign_gamma(double arg){
  double out;
  if(arg > 0){
    out = 1.0;
  }else{
    double integer = std::trunc(arg);
    out = std::pow(-1.0, integer + 1.0);
  }
  return(out);
}



Eigen::VectorXd lratio_gamma(double arg1, double arg2, double arg3, double arg4){
  Eigen::VectorXd out(2);
  double lratio_val = gsl_sf_lngamma(arg1) + gsl_sf_lngamma(arg2) -
    gsl_sf_lngamma(arg3) - gsl_sf_lngamma(arg4);
  double lratio_sign = sign_gamma(arg1)*sign_gamma(arg2)*sign_gamma(arg3)*sign_gamma(arg4);
  out[0] = lratio_val;
  out[1] = lratio_sign;

  return(out);
}




Eigen::VectorXd l1_e(double q_pri, double p_pri, double alpha,
                     double alpha_til, double s, double x){

  Eigen::VectorXd gs(4);

  f21_wrap(q_pri + alpha, q_pri + p_pri, q_pri + alpha_til,
           - 1.0 / (s * x), gs);
  if(gs[1]<0 | !std::isfinite(gs[1]) | gs[0] != 0) {
    gs[0] = 1.0;
    gs[1] = NAN;
    gs[2] = 1.0 / 0.0;
  } else {
    gs[1] = std::log(gs[1]);
    gs[2] = std::log(gs[2]) - gs[1];
  }

  return(gs);
}

Eigen::VectorXd l2_e(double q_pri, double p_pri, double alpha,
                     double alpha_til, double s, double x){
  Eigen::VectorXd gs(4);

  f21_wrap(alpha_til - alpha, alpha_til - p_pri, q_pri + alpha_til,
           - 1.0 / (s * x), gs);
  if(gs[1]<0 | !std::isfinite(gs[1]) | gs[0] != 0) {
    gs[0] = 1.0;
    gs[1] = NAN;
    gs[2] = 1.0 / 0.0;
  } else {
    gs[1] = std::log(gs[1]) +
      (alpha_til - alpha - p_pri - q_pri) * std::log(1.0 + 1.0/(s*x));
    gs[2] = std::log(gs[2]) +
      (alpha_til - alpha - p_pri - q_pri) * std::log(1.0 + 1.0/(s*x));
    gs[2] = gs[2] - gs[1];
  }

  return(gs);
}



Eigen::VectorXd l3_e(double q_pri, double p_pri, double alpha,
                     double alpha_til, double s, double x){
  Eigen::VectorXd gs(4);
  f21_wrap(q_pri + alpha, alpha_til - p_pri, q_pri + alpha_til,
           1.0 / (1.0 + s * x),gs);
  if(gs[1]<0 | !std::isfinite(gs[1]) | gs[0] != 0) {
    gs[0] = 1.0;
    gs[1] = NAN;
    gs[2] = 1.0 / 0.0;
  } else {
    gs[1] = std::log(gs[1]) +
      (- alpha - q_pri) * std::log(1.0 + 1.0/(s*x));
    gs[2] = std::log(gs[2]) +
      (- alpha - q_pri) * std::log(1.0 + 1.0/(s*x));
    gs[2] = gs[2] - gs[1];
  }

  return(gs);
}


Eigen::VectorXd l4_e(double q_pri, double p_pri, double alpha,
                     double alpha_til, double s, double x){

  Eigen::VectorXd gs(4);
  f21_wrap(q_pri + p_pri, alpha_til - alpha, q_pri + alpha_til,
           1.0 / (1.0 + s * x), gs);
  if(gs[1]<0 | !std::isfinite(gs[1]) | gs[0] != 0) {
    gs[0] = 1.0;
    gs[1] = NAN;
    gs[2] = 1.0 / 0.0;
  } else {
    gs[1] = std::log(gs[1]) +
      (- p_pri - q_pri) * std::log(1.0 + 1.0/(s*x));
    gs[2] = std::log(gs[2]) +
      (- p_pri - q_pri) * std::log(1.0 + 1.0/(s*x));
    gs[2] = gs[2] - gs[1];
  }

  return(gs);
}


Eigen::VectorXd l51_e(double q_pri, double p_pri, double alpha,
                      double alpha_til, double s, double x){

  Eigen::VectorXd gs(4);
  f21_wrap(q_pri + alpha, 1.0 + alpha - alpha_til, alpha - p_pri + 1.0,
           - (s * x), gs);
  if(gs[1]<0 | !std::isfinite(gs[1]) | gs[0] != 0) {
    gs[0] = 1.0;
    gs[1] = NAN;
    gs[2] = 1.0 / 0.0;
  } else {
    Eigen::VectorXd lrat = lratio_gamma(q_pri + alpha_til, p_pri - alpha,
                                        p_pri + q_pri, alpha_til - alpha);
    gs[1] = std::log(gs[1]) + lrat[0] +
      (q_pri + alpha) * std::log(s * x);
    gs[2] = std::log(gs[2]) + lrat[0] +
      (q_pri + alpha) * std::log(s * x);
    gs[2] = gs[2] - gs[1];
    gs[3] = lrat[1];
  }

  return(gs);
}


Eigen::VectorXd l52_e(double q_pri, double p_pri, double alpha,
                      double alpha_til, double s, double x){

  Eigen::VectorXd gs(4);
  f21_wrap(q_pri + p_pri, 1.0 + p_pri - alpha_til, p_pri - alpha + 1.0,
           - (s * x), gs);
  if(!std::isfinite(gs[1]) | gs[0] != 0) {
    gs[0] = 1.0;
    gs[1] = NAN;
    gs[2] = 1.0 / 0.0;
  } else {
    Eigen::VectorXd lrat = lratio_gamma(q_pri + alpha_til, alpha - p_pri,
                                        q_pri + alpha, alpha_til - p_pri);
    gs[3] <- sign(gs[1]) * 1.0;
    gs[1] = std::log(gs[1]) + lrat[0] +
      (q_pri + p_pri) * std::log(s * x);
    gs[2] = std::log(gs[2]) + lrat[0] +
      (q_pri + p_pri) * std::log(s * x);
    gs[2] = gs[2] - gs[1];
  }

  return(gs);
}


Eigen::VectorXd l61_e(double q_pri, double p_pri, double alpha,
                      double alpha_til, double s, double x){
  Eigen::VectorXd gs(4);
  f21_wrap(q_pri + alpha, alpha_til - p_pri, alpha - p_pri + 1.0,
           (s * x) / (1.0 + s * x), gs);
  if(gs[1]<0 | !std::isfinite(gs[1]) | gs[0] != 0) {
    gs[0] = 1.0;
    gs[1] = NAN;
    gs[2] = 1.0 / 0.0;
  } else {
    Eigen::VectorXd lrat = lratio_gamma(q_pri + alpha_til, p_pri - alpha,
                                        p_pri + q_pri, alpha_til - alpha);
    gs[1] = std::log(gs[1]) + lrat[0] +
      (- q_pri - alpha) * std::log(1.0 + 1.0 / (s * x));
    gs[2] = std::log(gs[2]) + lrat[0] +
      (- q_pri - alpha) * std::log(1.0 + 1.0 / (s * x));
    gs[2] = gs[2] - gs[1];
    gs[3] = lrat[1];
  }

  return(gs);
}


Eigen::VectorXd l62_e(double q_pri, double p_pri, double alpha,
                      double alpha_til, double s, double x){
  Eigen::VectorXd gs(4);
  f21_wrap(q_pri + p_pri, alpha_til - alpha, p_pri - alpha + 1.0,
           (s * x) / (1.0 + s * x), gs);
  if(!std::isfinite(gs[1]) | gs[0] != 0) {
    gs[0] = 1.0;
    gs[1] = NAN;
    gs[2] = 1.0 / 0.0;
  } else {
    Eigen::VectorXd lrat = lratio_gamma(q_pri + alpha_til, alpha - p_pri,
                                        q_pri + alpha, alpha_til - p_pri);
    gs[3] <- sign(gs[1]) * 1.0;
    gs[1] = std::log(gs[1]) + lrat[0] +
      (- q_pri - p_pri) * std::log(1.0 + 1.0 / (s * x));
    gs[2] = std::log(gs[2]) + lrat[0] +
      (- q_pri - p_pri) * std::log(1.0 + 1.0 / (s * x));
    gs[2] = gs[2] - gs[1];
  }

  return(gs);
}


Eigen::VectorXd l72_e(double q_pri, double p_pri, double alpha,
                      double alpha_til, double s, double x){

  Eigen::VectorXd gs(4);
  f21_wrap(alpha_til - alpha, 1.0 - alpha - q_pri, p_pri - alpha + 1.0,
           - (s * x), gs);
  if(!std::isfinite(gs[1]) | gs[0] != 0) {
    gs[0] = 1.0;
    gs[1] = NAN;
    gs[2] = 1.0 / 0.0;
  } else {
    Eigen::VectorXd lrat = lratio_gamma(q_pri + alpha_til, alpha - p_pri,
                                        q_pri + alpha, alpha_til - p_pri);
    gs[3] <- sign(gs[1]) * 1.0;
    gs[1] = std::log(gs[1]) + lrat[0] +
      (q_pri + p_pri) * std::log(s*x) + (alpha_til - alpha - p_pri - q_pri) * std::log(1.0 + s * x);
    gs[2] = std::log(gs[2]) + lrat[0] +
      (q_pri + p_pri) * std::log(s*x) + (alpha_til - alpha - p_pri - q_pri) * std::log(1.0 + s * x);
    gs[2] = gs[2] - gs[1];
  }

  return(gs);
}


Eigen::VectorXd l_sum_1(double q_pri, double p_pri,
                        double alpha, double alpha_til,
                        double s, double x){
  // calcola la 61 (salva val, err, status, sign)
  Eigen::VectorXd  l61 = l61_e(q_pri, p_pri, alpha, alpha_til, s, x);
  Eigen::VectorXd  l_1_out;
  if(s > 1.0/x){
    l_1_out = l61;
  }else{
    Eigen::VectorXd  l51 = l51_e(q_pri, p_pri, alpha, alpha_til, s, x);
    if(l51[2] > l61[2]){
      l_1_out = l61;
    }else{
      l_1_out = l51;
    }
  }

  return(l_1_out);
}


Eigen::VectorXd l_sum_2(double q_pri, double p_pri,
                        double alpha, double alpha_til,
                        double s, double x){
  Eigen::VectorXd  l62 = l62_e(q_pri, p_pri, alpha, alpha_til, s, x);
  Eigen::VectorXd  l_2_out;
  if(s > 1.0 / x){
    l_2_out = l62;
  }else{
    Eigen::VectorXd l52 = l52_e(q_pri, p_pri, alpha, alpha_til, s, x);
    Eigen::VectorXd l72 = l72_e(q_pri, p_pri, alpha, alpha_til, s, x);
    if(l52[2] > l62[2]){
      l_2_out = l62;
    }else{
      l_2_out = l52;
    }
    if(l_2_out[2] > l72[2]){
      l_2_out = l72;
    }

  }
  return(l_2_out);
}



Eigen::VectorXd l_sum(double q_pri, double p_pri,
                      double alpha, double alpha_til,
                      double s, double x, double lmlt){
  Eigen::VectorXd out(4);
  out[3] = 1.0;
  Eigen::VectorXd  l_1 = l_sum_1(q_pri, p_pri,
                                 alpha, alpha_til,
                                 s, x);
  Eigen::VectorXd  l_2 = l_sum_2(q_pri, p_pri,
                                 alpha, alpha_til,
                                 s, x);
  if(std::isinf(l_1[2]) | std::isinf(l_2[2])){
    out[1] = NAN;
    out[2] = 1.0/0.0;
    out[0] = 1;
  }else{
    if(s < 1/x){
      out[1] = std::log(l_1[3] * std::exp(l_1[1] + lmlt) +
        l_2[3] * std::exp(l_2[1] + lmlt)) - lmlt;
      out[2] = std::log(std::exp(l_1[2] + l_1[1] + lmlt) +
        std::exp(l_2[2] + l_2[1] + lmlt)) - out[1] - lmlt;
      out[0] = 0;
    }else{
      out[1] = std::log(l_1[3] * std::exp(l_1[1]) +
        l_2[3] * std::exp(l_2[1]));
      out[2] = std::log(std::exp(l_1[2] + l_1[1]) +
        std::exp(l_2[2] + l_2[1])) - out[1];
      out[0] = 0;
    }
  }

  return(out);

}




//'
// [[Rcpp::export]]

Eigen::VectorXd l_2F1(double q_pri, double p_pri,
                      double alpha, double alpha_til,
                      double s, double x, double lmlt){
  Eigen::VectorXd l_out;
  Eigen::VectorXd l_3 = l3_e(q_pri, p_pri, alpha, alpha_til, s, x);
  Eigen::VectorXd l_4 = l4_e(q_pri, p_pri, alpha, alpha_til, s, x);
  Eigen::VectorXd l_s = l_sum(q_pri, p_pri, alpha, alpha_til, s, x, lmlt);
  if(l_3[2] > l_4[2]){
    l_out = l_4;
  }else{
    l_out = l_3;
  }
  if(l_out[2] > l_s[2]){
    l_out = l_s;
  }

  if(s > 1.0 / x){
    Eigen::VectorXd l_1 = l1_e(q_pri, p_pri, alpha, alpha_til, s, x);
    Eigen::VectorXd l_2 = l2_e(q_pri, p_pri, alpha, alpha_til, s, x);
    if(l_out[2] > l_1[2]){
      l_out = l_1;
    }
    if(l_out[2] > l_2[2]){
      l_out = l_2;
    }
  }
  return(l_out);
}




