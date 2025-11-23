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
    x_noise = sensitivity::tell(gensol, Y_noise),
    x_qoi = sensitivity::tell(gensol, Y_qoi)
  )
}
