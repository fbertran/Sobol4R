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
