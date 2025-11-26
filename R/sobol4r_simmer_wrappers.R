#' M/M/1 queue model wrapper for Sobol designs
#'
#' Evaluate a simple M/M/1 queue built with \pkg{simmer} for each row of a
#' Sobol design matrix. The quantity of interest is the mean time in system
#' across \code{nrep} independent replications.
#'
#' @param X Design matrix or data.frame with columns \code{lambda} (arrival
#'   rate) and \code{mu} (service rate).
#' @param horizon Simulation horizon.
#' @param warmup Warmup period; arrivals ending before this time are discarded
#'   from the summary statistic.
#' @param nrep Number of replications used to average the mean time in system.
#' @return Numeric vector of length \code{nrow(X)}.
#' @export
sobol4r_mm1_model <- function(X, horizon = 1000, warmup = 200, nrep = 20L) {
  if (rlang::is_installed("simmer")) {
    stop("The 'simmer' package is required for sobol4r_mm1_model().", call. = FALSE)
  }
  X <- as.data.frame(X)
  if (!all(c("lambda", "mu") %in% names(X))) {
    stop("Input design must contain columns 'lambda' and 'mu'.")
  }
  nrep <- as.integer(nrep)
  if (nrep < 1L) {
    stop("nrep must be a positive integer.")
  }
  simulate_once <- function(lambda, mu) {
    traj <- simmer::trajectory("customer")
    traj <- simmer::seize(traj, "server")
    traj <- simmer::timeout(traj, function() stats::rexp(1, rate = mu))
    traj <- simmer::release(traj, "server")
    env <- simmer::simmer("MM1_queue")
    env <- simmer::add_resource(env, "server", capacity = 1, queue_size = Inf)
    env <- simmer::add_generator(env, "cust", traj,
                                 function() stats::rexp(1, rate = lambda))
    env <- simmer::run(env, until = horizon)
    arr <- simmer::get_mon_arrivals(env)
    arr <- arr[arr$end_time > warmup, , drop = FALSE]
    if (nrow(arr) == 0L) {
      return(0)
    }
    ts <- arr$end_time - arr$start_time
    m <- mean(ts, na.rm = TRUE)
    if (!is.finite(m)) 0 else m
  }
  vapply(seq_len(nrow(X)), function(i) {
    lambda <- X$lambda[i]
    mu <- X$mu[i]
    vals <- replicate(nrep, simulate_once(lambda, mu))
    mean(vals)
  }, numeric(1L))
}

#' Two-step clinic model wrapper for Sobol designs
#'
#' Simulate a simple clinic with separate registration and examination stages
#' using \pkg{simmer}. The quantity of interest is the mean time in system over
#' \code{nrep} replications for each parameter set.
#'
#' @param X Design matrix or data.frame with columns \code{lambda} (arrival
#'   rate), \code{mu_reg} (registration service rate), and \code{mu_exam}
#'   (examination service rate).
#' @param cap_reg,cap_exam Integer capacities for the registration and
#'   examination resources.
#' @param horizon Simulation horizon.
#' @param warmup_prob Fraction of the horizon treated as warmup and discarded
#'   before computing the mean time in system.
#' @param nrep Number of replications used to average the mean time in system.
#' @return Numeric vector of length \code{nrow(X)}.
#' @export
sobol4r_clinic_model <- function(X, cap_reg = 2, cap_exam = 3, horizon = 2000,
                                 warmup_prob = 0.2, nrep = 10L) {
  if (rlang::is_installed("simmer")) {
    stop("The 'simmer' package is required for sobol4r_clinic_model().", call. = FALSE)
  }
  X <- as.data.frame(X)
  if (!all(c("lambda", "mu_reg", "mu_exam") %in% names(X))) {
    stop("Input design must contain columns 'lambda', 'mu_reg', and 'mu_exam'.")
  }
  nrep <- as.integer(nrep)
  if (nrep < 1L) {
    stop("nrep must be a positive integer.")
  }
  warmup <- horizon * warmup_prob
  simulate_once <- function(lambda, mu_reg, mu_exam) {
    traj <- simmer::trajectory("patient")
    traj <- simmer::seize(traj, "registration")
    traj <- simmer::timeout(traj, function() stats::rexp(1, rate = mu_reg))
    traj <- simmer::release(traj, "registration")
    traj <- simmer::seize(traj, "examination")
    traj <- simmer::timeout(traj, function() stats::rexp(1, rate = mu_exam))
    traj <- simmer::release(traj, "examination")
    env <- simmer::simmer("clinic")
    env <- simmer::add_resource(env, "registration", capacity = cap_reg, queue_size = Inf)
    env <- simmer::add_resource(env, "examination", capacity = cap_exam, queue_size = Inf)
    env <- simmer::add_generator(env, "patient", traj,
                                 function() stats::rexp(1, rate = lambda))
    env <- simmer::run(env, until = horizon)
    arr <- simmer::get_mon_arrivals(env)
    arr <- arr[arr$end_time > warmup, , drop = FALSE]
    if (nrow(arr) == 0L) {
      return(0)
    }
    ts <- arr$end_time - arr$start_time
    m <- mean(ts, na.rm = TRUE)
    if (!is.finite(m)) 0 else m
  }
  vapply(seq_len(nrow(X)), function(i) {
    lambda <- X$lambda[i]
    mu_reg <- X$mu_reg[i]
    mu_exam <- X$mu_exam[i]
    vals <- replicate(nrep, simulate_once(lambda, mu_reg, mu_exam))
    mean(vals)
  }, numeric(1L))
}