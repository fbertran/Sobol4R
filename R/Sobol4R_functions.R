# Sobol4R - core functions extracted from Sobol_RV examples
# This file is meant to be placed in R/ inside the Sobol4R package.
# It refactors the original scripts into reusable functions.
#
# Imports (to be declared in DESCRIPTION / NAMESPACE):
#   Imports: stats, sensitivity
#
# The goal is to:
# - provide deterministic Sobol g functions
# - provide stochastic variants with additive noise on the output
# - provide variants with covariate dependent noise
# - provide helpers for a simple stochastic process model
# - provide a small helper to run Sobol indices with an optional QoI wrapper

#' Sobol g function (Saltelli) for arbitrary dimension
#'
#' Generic implementation of the non monotonic Sobol g function.
#'
#' @param X Numeric matrix or data.frame of inputs in [0,1].
#'   Rows are observations, columns are factors.
#' @param a Numeric vector of non negative parameters (importance factors).
#'   Its length must be at least the number of columns used.
#' @param use_cols Integer vector giving which columns of X to use.
#'   Defaults to all columns.
#'
#' @return Numeric vector of length nrow(X) with the model output.
#' @export
sobol4r_g_function <- function(X,
                               a = c(0, 1, 4.5, 9, 99, 99, 99, 99),
                               use_cols = seq_len(ncol(as.data.frame(X)))) {
  X <- as.data.frame(X)
  if (length(a) < length(use_cols)) {
    stop("Length of 'a' must be at least the number of used columns")
  }
  n <- nrow(X)
  if (n == 0L) return(numeric(0L))
  y <- rep(1, n)
  for (j in use_cols) {
    aj <- a[j]
    y <- y * (abs(4 * X[[j]] - 2) + aj) / (1 + aj)
  }
  y
}

#' Sobol g function restricted to the first two inputs
#'
#' This reproduces the function called `sobol.fun2` in the script.
#'
#' @inheritParams sobol4r_g_function
#' @return Numeric vector of model outputs.
#' @export
sobol4r_g2 <- function(X,
                       a = c(0, 1, 4.5, 9, 99, 99, 99, 99)) {
  sobol4r_g_function(X, a = a, use_cols = 1:2)
}

#' Additive output noise: constant Gaussian noise
#'
#' This corresponds to `sobol.fun3` in the script (noise N(0,1)).
#'
#' @param X Numeric matrix or data.frame of inputs.
#' @param sigma Standard deviation of the Gaussian noise.
#' @inheritParams sobol4r_g_function
#'
#' @return Numeric vector with g2(X) plus Gaussian noise.
#' @export
sobol4r_g2_noise_const <- function(X,
                                   a = c(0, 1, 4.5, 9, 99, 99, 99, 99),
                                   sigma = 1) {
  mu <- sobol4r_g2(X, a = a)
  mu + stats::rnorm(n = length(mu), mean = 0, sd = sigma)
}

#' Quantity of interest wrapper for stochastic output
#'
#' Compute a quantity of interest (QoI) from repeated evaluations
#' of a stochastic model, row wise.
#'
#' @param model Function taking a design matrix X and returning
#'   a numeric vector of length nrow(X). The model can be stochastic.
#' @param X Numeric matrix or data.frame of inputs.
#' @param n_rep Number of Monte Carlo replications for each row of X.
#' @param qoi Function that takes a numeric vector and returns a scalar
#'   (for example `mean`, `stats::var`, `stats::quantile` with `prob` fixed).
#' @param ... Additional arguments passed to `model`.
#'
#' @return Numeric vector of length nrow(X) with QoI values.
#' @export
sobol4r_qoi <- function(model, X,
                        n_rep = 1000L,
                        qoi = mean,
                        ...) {
  X <- as.data.frame(X)
  n <- nrow(X)
  if (n == 0L) return(numeric(0L))
  n_rep <- as.integer(n_rep)
  if (n_rep < 1L) stop("'n_rep' must be at least 1")
  
  # Accumulate QoI row wise
  out <- numeric(n)
  for (i in seq_len(n)) {
    vals <- model(X[i, , drop = FALSE], ...)
    if (n_rep > 1L) {
      extra <- replicate(n_rep - 1L, model(X[i, , drop = FALSE], ...))
      vals <- c(vals, extra)
    }
    out[i] <- qoi(vals)
  }
  out
}

#' Sobol g function with constant noise and QoI equal to the mean
#'
#' This corresponds to `sobol.fun4` in the script.
#'
#' @inheritParams sobol4r_g2_noise_const
#' @param n_rep Number of replications per design point for the QoI.
#'
#' @return Numeric vector of QoI values.
#' @export
sobol4r_g2_noise_const_qoi_mean <- function(X,
                                            a = c(0, 1, 4.5, 9, 99, 99, 99, 99),
                                            sigma = 1,
                                            n_rep = 1000L) {
  sobol4r_qoi(model = sobol4r_g2_noise_const,
              X = X,
              n_rep = n_rep,
              qoi = mean,
              a = a,
              sigma = sigma)
}

#' Sobol g function with covariate dependent noise
#'
#' This corresponds to `sobol.fun5` and `sobol.fun7` in the script,
#' where the noise is Gaussian with mean equal to one input column.
#'
#' @param X Numeric matrix or data.frame of inputs.
#' @param idx_mean Integer, index of the column that carries the mean of the noise.
#' @param noise_sd Standard deviation of the noise.
#' @inheritParams sobol4r_g_function
#'
#' @return Numeric vector of outputs.
#' @export
sobol4r_g2_noise_covariate <- function(X,
                                       a = c(0, 1, 4.5, 9, 99, 99, 99, 99),
                                       idx_mean = 3L,
                                       noise_sd = 1) {
  X <- as.data.frame(X)
  mu <- sobol4r_g2(X, a = a)
  mu + stats::rnorm(nrow(X), mean = X[[idx_mean]], sd = noise_sd)
}

#' QoI wrapper for the covariate dependent noise G function
#'
#' This corresponds to `sobol.fun6` and `sobol.fun8` in the script.
#'
#' @inheritParams sobol4r_g2_noise_covariate
#' @param n_rep Number of replications per design point for the QoI.
#'
#' @return Numeric vector of QoI values.
#' @export
sobol4r_g2_noise_covariate_qoi_mean <- function(X,
                                                a = c(0, 1, 4.5, 9, 99, 99, 99, 99),
                                                idx_mean = 3L,
                                                noise_sd = 1,
                                                n_rep = 1000L) {
  sobol4r_qoi(model = sobol4r_g2_noise_covariate,
              X = X,
              n_rep = n_rep,
              qoi = mean,
              a = a,
              idx_mean = idx_mean,
              noise_sd = noise_sd)
}

#' One unit of the simple stochastic process
#'
#' This is the `one_unit` function from the script.
#'
#' @param lambda1 Rate of the first exponential time.
#' @param lambda2 Rate of the second exponential time.
#' @param lambda3 Rate of the third exponential time.
#' @param p1 First Bernoulli success probability.
#' @param p2 Second Bernoulli success probability.
#'
#' @return Numeric vector of length 3:
#'   - indicator of success (0 or 1)
#'   - time to first event
#'   - additional time if success (0 otherwise)
#' @export
sobol4r_one_unit <- function(lambda1, lambda2, lambda3, p1, p2) {
  if (stats::rbinom(1, 1, p1) * stats::rbinom(1, 1, p2)) {
    c(
      1,
      stats::rexp(1, lambda1),
      stats::rexp(1, lambda2) + stats::rexp(1, lambda3)
    )
  } else {
    c(
      0,
      stats::rexp(1, lambda1),
      0
    )
  }
}

#' Time to M successes for one individual
#'
#' This is the `process.fun_indiv` function from the script.
#'
#' @param X_indiv Numeric vector of length 5:
#'   (lambda1, lambda2, lambda3, p1, p2).
#' @param M Target number of successes.
#'
#' @return Time to reach M successes for that individual.
#'   The value has an attribute "success" with the final number
#'   of successes (should be >= M).
#' @export
sobol4r_process_indiv <- function(X_indiv, M = 50L) {
  lambda1 <- X_indiv[1]
  lambda2 <- X_indiv[2]
  lambda3 <- X_indiv[3]
  p1 <- X_indiv[4]
  p2 <- X_indiv[5]
  
  success <- 0L
  time_start <- 0
  times_success <- numeric(0)
  
  while (success < M) {
    res_unit <- sobol4r_one_unit(lambda1, lambda2, lambda3, p1, p2)
    success <- success + res_unit[1]
    time_start <- time_start + res_unit[2]
    times_success <- c(
      times_success,
      res_unit[1] * (time_start + res_unit[3])
    )
  }
  time_Msuccess <- max(times_success)
  attr(time_Msuccess, "success") <- success
  time_Msuccess
}

#' Time to M successes for a design of individuals
#'
#' This is the `process.fun1` function from the script.
#'
#' @param X Numeric matrix or data.frame. Each row is one individual
#'   and must contain 5 columns (lambda1, lambda2, lambda3, p1, p2).
#' @param M Target number of successes.
#'
#' @return Numeric vector of length nrow(X) with the time to M successes.
#' @export
sobol4r_process <- function(X, M = 50L) {
  X <- as.matrix(X)
  apply(X, 1L, sobol4r_process_indiv, M = M)
}

#' QoI wrapper for the process model
#'
#' This is the `process.fun2` function from the script.
#'
#' @inheritParams sobol4r_process
#' @param n_rep Number of replications per individual for the QoI.
#' @param qoi Function that takes a numeric vector and returns a scalar.
#'
#' @return Numeric vector of length nrow(X) with QoI values.
#' @export
sobol4r_process_qoi <- function(X,
                                M = 50L,
                                n_rep = 10L,
                                qoi = mean) {
  X <- as.matrix(X)
  n <- nrow(X)
  if (n == 0L) return(numeric(0L))
  n_rep <- as.integer(n_rep)
  if (n_rep < 1L) stop("'n_rep' must be at least 1")
  
  out <- numeric(n)
  for (i in seq_len(n)) {
    vals <- sobol4r_process(X[i, , drop = FALSE], M = M)
    if (n_rep > 1L) {
      extra <- replicate(n_rep - 1L, sobol4r_process(X[i, , drop = FALSE], M = M))
      vals <- c(vals, extra)
    }
    out[i] <- qoi(vals)
  }
  out
}

#' Run Sobol analysis with optional QoI wrapper
#'
#' Helper around `sensitivity::sobol` that mimics the structure
#' of the original scripts. It never writes to disk.
#'
#' @param model Deterministic or stochastic model that takes a design X
#'   and returns a numeric vector of length nrow(X).
#' @param X1,X2 Matrices or data.frames used to build the Sobol design.
#' @param order Order of the Sobol indices (1 or 2).
#' @param nboot Number of bootstrap replicates for confidence intervals.
#' @param qoi Optional quantity of interest function. If not NULL,
#'   the model is evaluated repeatedly and QoI is computed row wise.
#' @param n_rep Number of replications per design row when `qoi` is not NULL.
#' @param ... Extra arguments passed to `model`.
#'
#' @return A `sobol` object (output of `sensitivity::tell`).
#' @export
sobol4r_run <- function(model,
                        X1,
                        X2,
                        order = 2,
                        nboot = 100L,
                        qoi = NULL,
                        n_rep = 1L,
                        ...) {
  if (!requireNamespace("sensitivity", quietly = TRUE)) {
    stop("Package 'sensitivity' is required but not installed")
  }
  sob_obj <- sensitivity::sobol(
    model = NULL,
    X1 = as.data.frame(X1),
    X2 = as.data.frame(X2),
    order = order,
    nboot = nboot
  )
  
  if (is.null(qoi) || n_rep <= 1L) {
    Y <- model(sob_obj$X, ...)
  } else {
    Y <- sobol4r_qoi(model = model,
                     X = sob_obj$X,
                     n_rep = n_rep,
                     qoi = qoi,
                     ...)
  }
  sensitivity::tell(sob_obj, Y)
}
















































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
sobol_g_function <- function(X,
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
sobol_g2_function <- function(X,
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
sobol_g2_additive_noise <- function(X, sd = 1) {
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
sobol_g2_qoi_mean <- function(X, nrep = 1000, sd = 1) {
  X <- as.matrix(X)
  # replicate returns a matrix of dimension (nrow(X) x nrep)
  sims <- replicate(nrep, sobol_g2_additive_noise(X, sd = sd))
  # rowMeans works whether sims is matrix or higher-dim array
  rowMeans(sims)
}

#' Design generation for Sobol indices
#'
#' Simple helper that wraps \code{sensitivity::sobol} with \code{model = NULL}
#' to create the extended design matrix used to evaluate the model.
#'
#' @param X1 First sample (matrix or data.frame).
#' @param X2 Second sample (matrix or data.frame).
#' @param order Maximum interaction order (1 or 2).
#' @param nboot Number of bootstrap replicates for confidence intervals.
#' @param ... Additional arguments passed to \code{sensitivity::sobol}.
#'
#' @return An object of class \code{"sobol"} whose \code{$X} field contains
#'   the design matrix. You should evaluate your model on \code{$X} and
#'   then call \code{sensitivity::tell()}.
#' @export
sobol4r_design <- function(X1, X2, order = 2, nboot = 0, ...) {
  sensitivity::sobol(
    model = NULL,
    X1 = X1,
    X2 = X2,
    order = order,
    nboot = nboot,
    ...
  )
}

#' Generic QoI-based Sobol indices for a stochastic model
#'
#' This function extends the classical Sobol indices to a stochastic
#' simulator by first computing a quantity of interest (QoI) for each
#' input point, such as the mean of repeated runs.
#'
#' @param model Stochastic model function that takes a matrix or data.frame
#'   \code{X} and returns a numeric vector of length \code{nrow(X)}.
#' @param X1,X2 Two base designs (matrices or data.frames).
#' @param qoi_fun Function used to summarize the repetitions
#'   (default is \code{mean}).
#' @param nrep Number of repetitions of the stochastic model for each
#'   design point.
#' @param order Maximum interaction order (1 or 2).
#' @param nboot Number of bootstrap replicates for Sobol indices.
#' @param ... Additional arguments passed to \code{model}.
#'
#' @return An object of class \code{"sobol"} with QoI-based Sobol indices.
#' @export
sobol_qoi_indices <- function(model,
                              X1,
                              X2,
                              qoi_fun = base::mean,
                              nrep = 1000,
                              order = 2,
                              nboot = 0,
                              ...) {
  gensol <- sobol4r_design(X1 = X1, X2 = X2, order = order, nboot = nboot)
  X_all <- gensol$X
  X_all <- as.matrix(X_all)
  
  # Evaluate the stochastic model nrep times
  sims <- replicate(nrep, model(X_all, ...))
  # sims has dimension (nrow(X_all) x nrep)
  if (!is.matrix(sims)) {
    sims <- as.matrix(sims)
  }
  
  # Apply QoI function row-wise
  qoi_vals <- apply(sims, 1L, qoi_fun)
  
  # Feed back into sensitivity::tell
  sensitivity::tell(gensol, qoi_vals)
}

# =============================================================================
# Example scenarios, following the Sobol4R slides and Rmd
# =============================================================================

#' Example 1: Deterministic G-function (reference case)
#'
#' Reproduces the classical non-random Sobol analysis on the G-function
#' with k = 8 inputs on \[0, 1\].
#'
#' @param n Monte Carlo sample size for each base design.
#' @param a Parameter vector for the G-function.
#' @param order Maximum interaction order for Sobol indices.
#' @param nboot Number of bootstrap replicates.
#'
#' @return An object of class \code{"sobol"}.
#' @export
sobol_example_g_deterministic <- function(n = 50000,
                                          a = c(0, 1, 4.5, 9, 99, 99, 99, 99),
                                          order = 2,
                                          nboot = 100) {
  X1 <- data.frame(matrix(stats::runif(8 * n), nrow = n))
  X2 <- data.frame(matrix(stats::runif(8 * n), nrow = n))
  
  gensol <- sobol4r_design(X1 = X1, X2 = X2, order = order, nboot = nboot)
  Y <- sobol_g_function(gensol$X, a = a)
  sensitivity::tell(gensol, Y)
}

#' Example 2: Random effect on the output (constant Gaussian noise)
#'
#' Two inputs in \[0, 1\], Sobol G-function with k = 2, plus additive
#' Gaussian noise, and a QoI based on the mean of repeated evaluations.
#'
#' @param n Monte Carlo sample size for each base design.
#' @param sd Standard deviation of the Gaussian noise.
#' @param nrep_qoi Number of repetitions for the QoI.
#' @param order Maximum interaction order.
#' @param nboot Number of bootstrap replicates.
#'
#' @return A list with three \code{"sobol"} objects:
#'   \code{x_det} (deterministic G-function),
#'   \code{x_noise} (single noisy output),
#'   \code{x_qoi} (QoI-based indices).
#' @export
sobol_example_random_output <- function(n = 50000,
                                        sd = 1,
                                        nrep_qoi = 1000,
                                        order = 2,
                                        nboot = 100) {
  X1 <- data.frame(
    C1 = stats::runif(n),
    C2 = stats::runif(n)
  )
  X2 <- data.frame(
    C1 = stats::runif(n),
    C2 = stats::runif(n)
  )
  
  gensol <- sobol4r_design(X1 = X1, X2 = X2, order = order, nboot = nboot)
  X_all <- gensol$X
  
  # Deterministic G-function with k = 2
  Y_det <- sobol_g2_function(X_all)
  
  # Single noisy realization
  Y_noise <- sobol_g2_additive_noise(X_all, sd = sd)
  
  # QoI: mean of nrep_qoi runs
  Y_qoi <- sobol_g2_qoi_mean(X_all, nrep = nrep_qoi, sd = sd)
  
  list(
    x_det = sensitivity::tell(gensol, Y_det),
    x_noise = sensitivity::tell(sobol4r_design(X1, X2, order, nboot), Y_noise),
    x_qoi = sensitivity::tell(sobol4r_design(X1, X2, order, nboot), Y_qoi)
  )
}

# Helpers for covariate dependent noise cases
sobol_g2_with_covariate_noise <- function(X) {
  # X with at least 3 columns, C1, C2, C3
  base <- sobol_g2_function(X)
  mu <- X[, 3]
  base + stats::rnorm(nrow(as.matrix(X)), mean = mu)
}

sobol_g2_qoi_covariate_mean <- function(X, nrep = 1000) {
  sims <- replicate(nrep, sobol_g2_with_covariate_noise(X))
  rowMeans(sims)
}

#' Example 3: Large covariate dependent random effect
#'
#' Third input C3 is uniform on \[1, 100\], used as the mean of a Gaussian
#' noise term added to the G-function. Quantity of interest is the mean
#' of repeated evaluations.
#'
#' @param n Monte Carlo sample size for each base design.
#' @param nrep_qoi Number of repetitions for the QoI.
#' @param order Maximum interaction order.
#' @param nboot Number of bootstrap replicates.
#'
#' @return A list with two \code{"sobol"} objects:
#'   \code{x_single} (single noisy run),
#'   \code{x_qoi} (QoI-based indices).
#' @export
sobol_example_covariate_large <- function(n = 50000,
                                          nrep_qoi = 1000,
                                          order = 2,
                                          nboot = 100) {
  X1 <- data.frame(
    C1 = stats::runif(n),
    C2 = stats::runif(n),
    C3 = stats::runif(n, min = 1, max = 100)
  )
  X2 <- data.frame(
    C1 = stats::runif(n),
    C2 = stats::runif(n),
    C3 = stats::runif(n, min = 1, max = 100)
  )
  
  gensol <- sobol4r_design(X1 = X1, X2 = X2, order = order, nboot = nboot)
  X_all <- gensol$X
  
  Y_single <- sobol_g2_with_covariate_noise(X_all)
  Y_qoi <- sobol_g2_qoi_covariate_mean(X_all, nrep = nrep_qoi)
  
  list(
    x_single = sensitivity::tell(gensol, Y_single),
    x_qoi = sensitivity::tell(sobol4r_design(X1, X2, order, nboot), Y_qoi)
  )
}

#' Example 4: Slight covariate dependent random effect
#'
#' Same as \code{sobol_example_covariate_large} but with C3 uniform
#' on \[1, 1.5\], that is with a much smaller range for the mean of the
#' Gaussian noise.
#'
#' @inheritParams sobol_example_covariate_large
#'
#' @return A list with two \code{"sobol"} objects:
#'   \code{x_single} (single noisy run),
#'   \code{x_qoi} (QoI-based indices).
#' @export
sobol_example_covariate_small <- function(n = 50000,
                                          nrep_qoi = 1000,
                                          order = 2,
                                          nboot = 100) {
  X1 <- data.frame(
    C1 = stats::runif(n),
    C2 = stats::runif(n),
    C3 = stats::runif(n, min = 1, max = 1.5)
  )
  X2 <- data.frame(
    C1 = stats::runif(n),
    C2 = stats::runif(n),
    C3 = stats::runif(n, min = 1, max = 1.5)
  )
  
  gensol <- sobol4r_design(X1 = X1, X2 = X2, order = order, nboot = nboot)
  X_all <- gensol$X
  
  Y_single <- sobol_g2_with_covariate_noise(X_all)
  Y_qoi <- sobol_g2_qoi_covariate_mean(X_all, nrep = nrep_qoi)
  
  list(
    x_single = sensitivity::tell(gensol, Y_single),
    x_qoi = sensitivity::tell(sobol4r_design(X1, X2, order, nboot), Y_qoi)
  )
}

# =============================================================================
# Simple process example (random distributional parameters)
# =============================================================================

#' Simulate one unit in the simple process
#'
#' Internal helper. The output is a length-3 numeric vector:
#' success indicator, waiting time to decision, and extra time if success.
#'
#' @param lambda1,lambda2,lambda3 Positive rates.
#' @param p1,p2 Success probabilities.
#'
#' @return Numeric vector c(success, t_decision, extra_time_if_success).
#' @keywords internal
one_unit <- function(lambda1, lambda2, lambda3, p1, p2) {
  success_flag <- stats::rbinom(1L, 1L, p1) * stats::rbinom(1L, 1L, p2)
  if (success_flag == 1L) {
    c(
      1,
      stats::rexp(1L, lambda1),
      stats::rexp(1L, lambda2) + stats::rexp(1L, lambda3)
    )
  } else {
    c(
      0,
      stats::rexp(1L, lambda1),
      0
    )
  }
}

#' Time to M successes for one individual
#'
#' Stochastic model that simulates successive units until M successes
#' occur, and returns the time when the M-th success happens.
#'
#' @param X_indiv Numeric vector \code{c(lambda1, lambda2, lambda3, p1, p2)}.
#' @param M Target number of successes.
#'
#' @return Scalar time to M successes, with attribute \code{"success"}.
#' @export
process_fun_indiv <- function(X_indiv, M = 50) {
  lambda1 <- X_indiv[1]
  lambda2 <- X_indiv[2]
  lambda3 <- X_indiv[3]
  p1 <- X_indiv[4]
  p2 <- X_indiv[5]
  
  success <- 0
  time_start <- 0
  times_success <- numeric(0)
  
  while (success < M) {
    res_unit <- one_unit(lambda1, lambda2, lambda3, p1, p2)
    success <- success + res_unit[1]
    time_start <- time_start + res_unit[2]
    times_success <- c(
      times_success,
      res_unit[1] * (time_start + res_unit[3])
    )
  }
  
  time_Msuccess <- max(times_success)
  attr(time_Msuccess, "success") <- success
  time_Msuccess
}

#' Process model for a matrix of individuals
#'
#' Applies \code{process_fun_indiv} row-wise to a matrix of parameters.
#'
#' @param X Matrix or data.frame with columns
#'   \code{lambda1, lambda2, lambda3, p1, p2}.
#' @param M Target number of successes.
#'
#' @return Numeric vector of length \code{nrow(X)}.
#' @export
process_fun1 <- function(X, M = 50) {
  X <- as.matrix(X)
  apply(X, 1L, process_fun_indiv, M = M)
}

#' QoI wrapper for the process model
#'
#' For each row of \code{X}, evaluates \code{process_fun1} several times
#' and returns the mean time to M successes.
#'
#' @param X Matrix or data.frame of parameters.
#' @param M Target number of successes.
#' @param nrep Number of repetitions for the QoI.
#'
#' @return Numeric vector of QoI values.
#' @export
process_fun2 <- function(X, M = 50, nrep = 10) {
  X <- as.matrix(X)
  sims <- replicate(nrep, process_fun1(X = X, M = M))
  rowMeans(sims)
}

#' Example 5: Sobol indices for the process model
#'
#' Computes Sobol indices for the simple process example with random
#' distributional parameters. Uses both a single trajectory and a
#' QoI based on repeated runs.
#'
#' @param n Monte Carlo sample size for each base design.
#' @param M Target number of successes.
#' @param nrep_qoi Number of repetitions for the QoI.
#' @param order Maximum interaction order.
#' @param nboot Number of bootstrap replicates.
#'
#' @return A list with two \code{"sobol"} objects:
#'   \code{xp_single} and \code{xp_qoi}.
#' @export
sobol_example_process <- function(n = 100,
                                  M = 50,
                                  nrep_qoi = 10,
                                  order = 1,
                                  nboot = 10) {
  # Draw random distributional parameters for each individual
  X1 <- data.frame(t(replicate(
    n,
    c(
      1 / stats::runif(1, min = 20, max = 100),
      1 / stats::runif(1, min = 24, max = 2000),
      1 / stats::runif(1, min = 24, max = 120),
      stats::runif(1, min = 0.05, max = 0.3),
      stats::runif(1, min = 0.3, max = 0.7)
    )
  )))
  X2 <- data.frame(t(replicate(
    n,
    c(
      1 / stats::runif(1, min = 20, max = 100),
      1 / stats::runif(1, min = 24, max = 2000),
      1 / stats::runif(1, min = 24, max = 120),
      stats::runif(1, min = 0.05, max = 0.3),
      stats::runif(1, min = 0.3, max = 0.7)
    )
  )))
  
  gensol <- sobol4r_design(X1 = X1, X2 = X2, order = order, nboot = nboot)
  X_all <- gensol$X
  
  Y_single <- process_fun1(X_all, M = M)
  Y_qoi <- process_fun2(X_all, M = M, nrep = nrep_qoi)
  
  list(
    xp_single = sensitivity::tell(gensol, Y_single),
    xp_qoi = sensitivity::tell(sobol4r_design(X1, X2, order, nboot), Y_qoi)
  )
}