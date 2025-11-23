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
  one_unit_cpp(lambda1, lambda2, lambda3, p1, p2)
}

one_unit_R <- function(lambda1, lambda2, lambda3, p1, p2) {
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
  process_fun_indiv_cpp(X_indiv, M)
}

process_fun_indiv_R <- function(X_indiv, M = 50) {
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
process_fun_row_wise <- function(X, M = 50) {
  process_fun_row_wise_cpp(as.matrix(X), M)
}

process_fun_row_wise_R <- function(X, M = 50) {
  X <- as.matrix(X)
  apply(X = X, MARGIN=1L, FUN=process_fun_indiv, M = M)
}

#' QoI wrapper for the process model
#'
#' For each row of \code{X}, evaluates \code{process_fun_row_wise} several times
#' and returns the mean time to M successes.
#'
#' @param X Matrix or data.frame of parameters.
#' @param M Target number of successes.
#' @param nrep Number of repetitions for the QoI.
#'
#' @return Numeric vector of QoI values.
#' @export
process_fun_mean_to_M <- function(X, M = 50, nrep = 10) {
  process_fun_mean_to_M_cpp(as.matrix(X), M, nrep)
}

process_fun_mean_to_M_R <- function(X, M = 50, nrep = 10) {
  X <- as.matrix(X)
  sims <- replicate(nrep, process_fun_row_wise(X = X, M = M))
  rowMeans(sims)
}