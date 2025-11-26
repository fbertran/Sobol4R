#' Sobol G-function (Saltelli reference function) - C++ backend
#'
#' Generic implementation of the Sobol G-function for k inputs.
#' Columns of \code{X} are interpreted as inputs X1, X2, ..., Xk.
#'
#' @param X Numeric matrix or data.frame of inputs in \[0, 1\].
#' @param a Numeric vector of parameters a_j controlling importance.
#'   Its length must be at least \code{ncol(X)}.
#'
#' @return Numeric vector of length \code{nrow(X)} with model outputs.
#' @export
sobol_g_function <- function(X,
                             a = c(0, 1, 4.5, 9, 99, 99, 99, 99)) {
  X <- as.matrix(X)
  sobol_g_cpp(X, a)
}

#' Sobol G-function restricted to the first two inputs - C++ backend
#' 
#' Convenience wrapper around \code{sobol_g_function} that uses
#' only the first two columns of \code{X}.
#'
#' @param X Numeric matrix or data.frame with at least two columns.
#' @param a Numeric vector of parameters (at least length 2).
#'
#' @return Numeric vector of length \code{nrow(X)} with model outputs.
#'
#' @export
sobol_g2_function <- function(X,
                              a = c(0, 1, 4.5, 9, 99, 99, 99, 99)) {
  X <- as.matrix(X)
  if (ncol(X) < 2L) stop("X must have at least two columns.")
  sobol_g2_cpp(X, a)
}

#' Additive Gaussian noise on the Sobol G-function (k = 2) - C++ backend
#' 
#' @param X Numeric matrix or data.frame with at least two columns.
#' @param sd Standard deviation of the Gaussian noise.
#' @param a Numeric vector of parameters (at least length 2).
#'
#' @return Numeric vector of model outputs with noise.
#' 
#' @export
sobol_g2_additive_noise <- function(X, sd = 1,
                                    a = c(0, 1, 4.5, 9, 99, 99, 99, 99)) {
  X <- as.matrix(X)
  sobol_g2_additive_noise_cpp(X, sd = sd, a = a)
}


#' QoI wrapper for the noisy G-function (k = 2) - C++ backend
#'
#' Computes a mean over repeated evaluations of the noisy model.
#'
#' @param X Numeric matrix or data.frame with at least two columns.
#' @param nrep Number of replicates used for the QoI.
#' @param sd Standard deviation of the Gaussian noise.
#' @param a Numeric vector of parameters (at least length 2).
#'
#' @return Numeric vector of QoI values (means over \code{nrep} runs).
#' @export
sobol_g2_qoi_mean <- function(X, nrep = 1000, sd = 1,
                              a = c(0, 1, 4.5, 9, 99, 99, 99, 99)) {
  X <- as.matrix(X)
  sobol_g2_qoi_mean_cpp(X, nrep = nrep, sd = sd, a = a)
}

#' Covariate dependent Gaussian noise on the Sobol G-function (k = 2) - C++ backend
#'
#' @param X Numeric matrix or data.frame with at least two columns.
#' @param a Numeric vector of parameters (at least length 2).
#'
#' @return Numeric vector of model outputs with noise.
#' 
#' @export
sobol_g2_with_covariate_noise <- function(X,
                                          a = c(0, 1, 4.5, 9, 99, 99, 99, 99)) {
  X <- as.matrix(X)
  sobol_g2_with_covariate_noise_cpp(X, a = a)
}

#' QoI wrapper for covariate noisy G-function (k = 2) - C++ backend
#'
#' Computes a mean over repeated evaluations of the noisy model.
#'
#' @param X Numeric matrix or data.frame with at least two columns.
#' @param nrep Number of replicates used for the QoI.
#' @param a Numeric vector of parameters (at least length 2).
#'
#' @return Numeric vector of QoI values (means over \code{nrep} runs).
#' 
#' @export
sobol_g2_qoi_covariate_mean <- function(X, nrep = 1000,
                                        a = c(0, 1, 4.5, 9, 99, 99, 99, 99)) {
  X <- as.matrix(X)
  sobol_g2_qoi_covariate_mean_cpp(X, nrep = nrep, a = a)
}


#' Fast Ishigami Test Function
#'
#' C++ implementation of the Ishigami function that is widely used as a
#' benchmark for Sobol sensitivity indices. The implementation is vectorised and
#' therefore convenient for Monte Carlo experiments.
#'
#' @param x Numeric matrix with three columns representing the inputs.
#' @param a Numeric scalar controlling the nonlinear term.
#' @param b Numeric scalar controlling the interaction term.
#' @return Numeric vector of simulator outputs.
#' @export
#' @examples
#' x <- matrix(runif(30, -pi, pi), ncol = 3)
#' ishigami_model(x)
ishigami_model <- function(x, a = 7, b = 0.1) {
  if (!is.matrix(x) || ncol(x) != 3) {
    stop("x must be a numeric matrix with three columns.")
  }
  ishigami_cpp( x[, 1], x[, 2], x[, 3], as.numeric(a), as.numeric(b))
}

