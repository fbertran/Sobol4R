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
    x_qoi = sensitivity::tell(gensol, Y_qoi)
  )
}