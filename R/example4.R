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
    x_qoi = sensitivity::tell(gensol, Y_qoi)
  )
}