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
  
  Y_single <- process_fun_row_wise(X_all, M = M)
  Y_qoi <- process_fun_mean_to_M(X_all, M = M, nrep = nrep_qoi)
  
  list(
    xp_single = sensitivity::tell(gensol, Y_single),
    xp_qoi = sensitivity::tell(gensol, Y_qoi)
  )
}