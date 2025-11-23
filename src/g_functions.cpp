// Fast C++ implementations of Sobol G-functions and noisy variants

#include <Rcpp.h>
using namespace Rcpp;

// Core Sobol G-function ----------------------------------------------

// [[Rcpp::export]]
NumericVector sobol_g_cpp(NumericMatrix X, NumericVector a) {
  const int n = X.nrow();
  const int k = X.ncol();
  
  if (a.size() < k) {
    stop("Length of 'a' must be at least ncol(X).");
  }
  
  NumericVector y(n, 1.0);
  
  for (int j = 0; j < k; ++j) {
    const double aj = a[j];
    const double denom = 1.0 + aj;
    for (int i = 0; i < n; ++i) {
      const double xij = X(i, j);
      const double term = (std::fabs(4.0 * xij - 2.0) + aj) / denom;
      y[i] *= term;
    }
  }
  return y;
}

// G-function restricted to first two inputs --------------------------

// [[Rcpp::export]]
NumericVector sobol_g2_cpp(NumericMatrix X, NumericVector a) {
  const int n = X.nrow();
  const int k = X.ncol();
  
  if (k < 2) {
    stop("X must have at least two columns.");
  }
  if (a.size() < 2) {
    stop("Length of 'a' must be at least 2.");
  }
  
  NumericVector y(n, 1.0);
  for (int j = 0; j < 2; ++j) {
    const double aj = a[j];
    const double denom = 1.0 + aj;
    for (int i = 0; i < n; ++i) {
      const double xij = X(i, j);
      const double term = (std::fabs(4.0 * xij - 2.0) + aj) / denom;
      y[i] *= term;
    }
  }
  return y;
}

// Additive Gaussian noise on G2 --------------------------------------

// [[Rcpp::export]]
NumericVector sobol_g2_additive_noise_cpp(NumericMatrix X,
                                          double sd = 1.0,
                                          NumericVector a = NumericVector::create(0.0, 1.0, 4.5, 9.0, 99.0, 99.0, 99.0, 99.0)) {
  const int n = X.nrow();
  if (X.ncol() < 2) {
    stop("X must have at least two columns.");
  }
  if (a.size() < 2) {
    stop("Length of 'a' must be at least 2.");
  }
  
  NumericVector base = sobol_g2_cpp(X, a);
  NumericVector out(n);
  
  for (int i = 0; i < n; ++i) {
    // mean 0, sd sd
    double noise = R::rnorm(0.0, sd);
    out[i] = base[i] + noise;
  }
  return out;
}

// QoI wrapper: mean of noisy G2 over nrep replicates -----------------

// [[Rcpp::export]]
NumericVector sobol_g2_qoi_mean_cpp(NumericMatrix X,
                                    int nrep = 1000,
                                    double sd = 1.0,
                                    NumericVector a = NumericVector::create(0.0, 1.0, 4.5, 9.0, 99.0, 99.0, 99.0, 99.0)) {
  const int n = X.nrow();
  if (X.ncol() < 2) {
    stop("X must have at least two columns.");
  }
  if (a.size() < 2) {
    stop("Length of 'a' must be at least 2.");
  }
  if (nrep < 1) {
    stop("nrep must be at least 1.");
  }
  
  // Deterministic part of the model
  NumericVector base = sobol_g2_cpp(X, a);
  NumericVector out(n);
  
  for (int i = 0; i < n; ++i) {
    double sum_noise = 0.0;
    for (int r = 0; r < nrep; ++r) {
      sum_noise += R::rnorm(0.0, sd);
    }
    double mean_noise = sum_noise / static_cast<double>(nrep);
    out[i] = base[i] + mean_noise;
  }
  return out;
}

// Covariate dependent noise on G2 ------------------------------------

// [[Rcpp::export]]
NumericVector sobol_g2_with_covariate_noise_cpp(NumericMatrix X,
                                                NumericVector a = NumericVector::create(0.0, 1.0, 4.5, 9.0, 99.0, 99.0, 99.0, 99.0)) {
  const int n = X.nrow();
  const int k = X.ncol();
  
  if (k < 3) {
    stop("X must have at least three columns for covariate noise.");
  }
  if (a.size() < 2) {
    stop("Length of 'a' must be at least 2.");
  }
  
  NumericVector base = sobol_g2_cpp(X, a);
  NumericVector out(n);
  
  for (int i = 0; i < n; ++i) {
    double mu = X(i, 2); // third column, zero based index 2
    double noise = R::rnorm(mu, 1.0);
    out[i] = base[i] + noise;
  }
  return out;
}

// QoI wrapper for covariate noise case -------------------------------

// [[Rcpp::export]]
NumericVector sobol_g2_qoi_covariate_mean_cpp(NumericMatrix X,
                                              int nrep = 1000,
                                              NumericVector a = NumericVector::create(0.0, 1.0, 4.5, 9.0, 99.0, 99.0, 99.0, 99.0)) {
  const int n = X.nrow();
  const int k = X.ncol();
  
  if (k < 3) {
    stop("X must have at least three columns for covariate noise.");
  }
  if (a.size() < 2) {
    stop("Length of 'a' must be at least 2.");
  }
  if (nrep < 1) {
    stop("nrep must be at least 1.");
  }
  
  NumericVector base = sobol_g2_cpp(X, a);
  NumericVector out(n);
  
  for (int i = 0; i < n; ++i) {
    double mu = X(i, 2); // third column
    double sum_noise = 0.0;
    for (int r = 0; r < nrep; ++r) {
      sum_noise += R::rnorm(mu, 1.0);
    }
    double mean_noise = sum_noise / static_cast<double>(nrep);
    out[i] = base[i] + mean_noise;
  }
  return out;
}