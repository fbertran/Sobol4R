# Sobol4R - core functions and examples
# This file is a proposed R/ script for the Sobol4R package.

#' Sobol G-function (Saltelli reference function)
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
sobol_g_function_R <- function(X,
                             a = c(0, 1, 4.5, 9, 99, 99, 99, 99)) {
  X <- as.matrix(X)
  k <- ncol(X)
  if (length(a) < k) {
    stop("Length of 'a' must be at least ncol(X).")
  }
  y <- rep(1, nrow(X))
  for (j in seq_len(k)) {
    y <- y * (abs(4 * X[, j] - 2) + a[j]) / (1 + a[j])
  }
  y
}

#' Sobol G-function restricted to the first two inputs
#'
#' Convenience wrapper around \code{sobol_g_function} that uses
#' only the first two columns of \code{X}.
#'
#' @param X Numeric matrix or data.frame with at least two columns.
#' @param a Numeric vector of parameters (at least length 2).
#'
#' @return Numeric vector of model outputs.
#' @export
sobol_g2_function_R <- function(X,
                              a = c(0, 1, 4.5, 9, 99, 99, 99, 99)) {
  X <- as.matrix(X)
  if (ncol(X) < 2) {
    stop("X must have at least two columns.")
  }
  sobol_g_function(X[, 1:2, drop = FALSE], a = a)
}

#' Additive Gaussian noise on the Sobol G-function (k = 2)
#'
#' @param X Numeric matrix or data.frame with at least two columns.
#' @param sd Standard deviation of the Gaussian noise.
#'
#' @return Numeric vector of model outputs with noise.
#' @export
sobol_g2_additive_noise_R <- function(X, sd = 1) {
  base <- sobol_g2_function(X)
  base + stats::rnorm(nrow(as.matrix(X)), sd = sd)
}

#' Quantity-of-interest wrapper for the noisy G-function (k = 2)
#'
#' Computes a mean over repeated evaluations of the noisy model.
#'
#' @param X Numeric matrix or data.frame with at least two columns.
#' @param nrep Number of replicates used for the QoI.
#' @param sd Standard deviation of the Gaussian noise.
#'
#' @return Numeric vector of QoI values (means over \code{nrep} runs).
#' @export
sobol_g2_qoi_mean_R <- function(X, nrep = 1000, sd = 1) {
  X <- as.matrix(X)
  # replicate returns a matrix of dimension (nrow(X) x nrep)
  sims <- replicate(nrep, sobol_g2_additive_noise(X, sd = sd))
  # rowMeans works whether sims is matrix or higher-dim array
  rowMeans(sims)
}

# Helpers for covariate dependent noise cases
#' Additive Gaussian noise on the Sobol G-function (k = 2)
#'
#' @param X Numeric matrix or data.frame with at least two columns.
#'
#' @return Numeric vector of model outputs with noise.
#' @export
sobol_g2_with_covariate_noise_R <- function(X) {
  # X with at least 3 columns, C1, C2, C3
  base <- sobol_g2_function(X)
  mu <- X[, 3]
  base + stats::rnorm(nrow(as.matrix(X)), mean = mu)
}

#' Quantity-of-interest wrapper for the covariate noisy G-function (k = 2)
#'
#' Computes a mean over repeated evaluations of the noisy model.
#'
#' @param X Numeric matrix or data.frame with at least two columns.
#' @param nrep Number of replicates used for the QoI.
#'
#' @return Numeric vector of QoI values (means over \code{nrep} runs).
#' @export
sobol_g2_qoi_covariate_mean_R <- function(X, nrep = 1000) {
  sims <- replicate(nrep, sobol_g2_with_covariate_noise(X))
  rowMeans(sims)
}