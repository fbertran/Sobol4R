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
                        qoi_fun = NULL,
                        nrep = 1L,
                        type = c("sobol","sobol2007"),
                        ...) {
  if (!requireNamespace("sensitivity", quietly = TRUE)) {
    stop("Package 'sensitivity' is required but not installed")
  }
  if (is.null(qoi_fun) || nrep <= 1L) {
    sob_obj <- sobol4r_design(
    X1 = as.data.frame(X1),
    X2 = as.data.frame(X2),
    order = order,
    nboot = nboot
  )
    Y <- model(sob_obj$X, ...)
    sensitivity::tell(sob_obj, Y)
  } else {
    Y <- sobol4r_qoi_indices(model = model,
                     X1 = as.data.frame(X1),
                     X2 = as.data.frame(X2),
                     qoi_fun = qoi_fun,
                     nrep = nrep,
                     order = order,
                     nboot = nboot,
                     type = type,
                     ...)
  }
  return(Y)
}
