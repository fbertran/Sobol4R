// Ishigami function in Rcpp
#include <Rcpp.h>
#include <cmath>

// [[Rcpp::export]]
Rcpp::NumericVector ishigami_cpp(Rcpp::NumericVector x1,
                                      Rcpp::NumericVector x2,
                                      Rcpp::NumericVector x3,
                                      double a,
                                      double b) {
  R_xlen_t n = x1.size();
  if (x2.size() != n || x3.size() != n) {
    Rcpp::stop("All inputs must share the same length");
  }
  
  Rcpp::NumericVector out(n);
  
  for (R_xlen_t i = 0; i < n; ++i) {
    double x1i = x1[i];
    double x2i = x2[i];
    double x3i = x3[i];
    
    double s1   = std::sin(x1i);
    double s2   = std::sin(x2i);
    double s2_2 = s2 * s2;
    double x3_2 = x3i * x3i;
    double x3_4 = x3_2 * x3_2;
    
    out[i] = s1 + a * s2_2 + b * x3_4 * s1;
  }
  
  return out;
}











