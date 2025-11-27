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
#' @param qoi_fun Optional quantity of interest function. If not NULL,
#'   the model is evaluated repeatedly and QoI is computed row wise.
#' @param nrep Number of replications per design row when `qoi` is not NULL.
#' @param type Type of Monte Carlo Estimation of Sobol' Indices to be used.
#'   Supported estimators mirror the \pkg{sensitivity} helpers: \code{sobol},
#'   \code{sobol2007}, \code{soboljansen}, \code{sobolEff}, and
#'   \code{sobolmartinez}. Defaults to \code{"soboljansen"} because it offers
#'   robust first and total order indices on both centred and non-centred
#'   outputs.
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
                        type = c("soboljansen", "sobol","sobol2007", "sobolEff", "sobolmartinez"),
                        ...) {
  type = match.arg(type, c("soboljansen", "sobol", "sobol2007", "sobolEff", "sobolmartinez"))
  if (!rlang::is_installed("sensitivity")) {
    stop("Package 'sensitivity' is required but not installed")
  }
  if (is.null(qoi_fun) || nrep <= 1L) {
    if(type == "sobol2007"){    
      sob_obj <- sobol4r_design(
      X1 = as.data.frame(X1),
      X2 = as.data.frame(X2),
      nboot = nboot,
      type = type)
      } else {
    sob_obj <- sobol4r_design(
    X1 = as.data.frame(X1),
    X2 = as.data.frame(X2),
    order = order,
    nboot = nboot,
    type = type)
    }
    Y_raw <- model(sob_obj$X, ...)
    Y <- Y_raw - mean(Y_raw)
    sensitivity::tell(sob_obj, Y)
  } else {
    sob_obj <- sobol4r_qoi_indices(model = model,
                     X1 = as.data.frame(X1),
                     X2 = as.data.frame(X2),
                     qoi_fun = qoi_fun,
                     nrep = nrep,
                     order = order,
                     nboot = nboot,
                     type = type,
                     ...)
  }
  return(sob_obj)
}
