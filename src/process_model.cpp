#include <Rcpp.h>
using namespace Rcpp;

/*
 Simulate one unit.
 Equivalent to R version:
 
 success_flag <- rbinom(1,1,p1)*rbinom(1,1,p2)
 if success: c(1, rexp(1,lambda1), rexp(1,lambda2)+rexp(1,lambda3))
 else:       c(0, rexp(1,lambda1), 0)
 */

// [[Rcpp::export]]
NumericVector one_unit_cpp(double lambda1,
                           double lambda2,
                           double lambda3,
                           double p1,
                           double p2) {
  int s1 = R::rbinom(1.0, p1);
  int s2 = R::rbinom(1.0, p2);
  int success_flag = s1 * s2;
  
  double t1 = R::rexp(lambda1);
  double t_extra = 0.0;
  
  if (success_flag == 1) {
    t_extra = R::rexp(lambda2) + R::rexp(lambda3);
  }
  
  NumericVector out(3);
  out[0] = success_flag;
  out[1] = t1;
  out[2] = t_extra;
  return out;
}

/*
 Simulate time to M successes for ONE individual (lambda1, lambda2, lambda3, p1, p2)
 */

// [[Rcpp::export]]
double process_fun_indiv_cpp(NumericVector x_indiv, int M = 50) {
  
  double lambda1 = x_indiv[0];
  double lambda2 = x_indiv[1];
  double lambda3 = x_indiv[2];
  double p1      = x_indiv[3];
  double p2      = x_indiv[4];
  
  int success = 0;
  double time_start = 0.0;
  double time_Msuccess = 0.0;
  
  while (success < M) {
    int s1 = R::rbinom(1.0, p1);
    int s2 = R::rbinom(1.0, p2);
    int success_flag = s1 * s2;
    
    double t1 = R::rexp(lambda1);
    time_start += t1;
    
    if (success_flag == 1) {
      double t_extra = R::rexp(lambda2) + R::rexp(lambda3);
      double t_success = time_start + t_extra;
      if (t_success > time_Msuccess) {
        time_Msuccess = t_success;
      }
      success++;
    }
  }
  
  return time_Msuccess;
}

/*
 Vectorised version: apply process_fun_indiv_cpp to each row of X.
 X is n by 5.
 */

// [[Rcpp::export]]
NumericVector process_fun_row_wise_cpp(NumericMatrix X, int M = 50) {
  int n = X.nrow();
  NumericVector out(n);
  
  for (int i = 0; i < n; i++) {
    NumericVector row = X(i, _);
    out[i] = process_fun_indiv_cpp(row, M);
  }
  return out;
}

/*
 QoI model: mean time to M successes over nrep simulations.
 */

// [[Rcpp::export]]
NumericVector process_fun_mean_to_M_cpp(NumericMatrix X,
                                        int M   = 50,
                                        int nrep = 10) {
  int n = X.nrow();
  NumericVector out(n);
  
  for (int i = 0; i < n; i++) {
    NumericVector row = X(i, _);
    double sum = 0.0;
    for (int r = 0; r < nrep; r++) {
      sum += process_fun_indiv_cpp(row, M);
    }
    out[i] = sum / nrep;
  }
  return out;
}
