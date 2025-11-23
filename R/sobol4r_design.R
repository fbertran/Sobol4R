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

