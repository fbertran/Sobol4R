#' Generic QoI-based Sobol indices for a stochastic model
#'
#' This function extends the classical Sobol indices to a stochastic
#' simulator by first computing a quantity of interest (QoI) for each
#' input point, such as the mean of repeated runs.
#'
#' @param model Stochastic model function that takes a matrix or data.frame
#'   \code{X} and returns a numeric vector of length \code{nrow(X)}.
#' @param X1,X2 Two base designs (matrices or data.frames).
#' @param qoi_fun Function used to summarize the repetitions
#'   (default is \code{mean}).
#' @param nrep Number of repetitions of the stochastic model for each
#'   design point.
#' @param order Maximum interaction order (1 or 2).
#' @param nboot Number of bootstrap replicates for Sobol indices.
#' @param ... Additional arguments passed to \code{model}.
#'
#' @return An object of class \code{"sobol"} with QoI-based Sobol indices.
#' @export
sobol4r_qoi_indices <- function(model,
                                X1,
                                X2,
                                qoi_fun = base::mean,
                                nrep    = 1000,
                                order   = 2,
                                nboot   = 0,
                                ...) {
  gensol <- sobol4r_design(X1 = X1, X2 = X2, order = order, nboot = nboot)
  X_all  <- as.matrix(gensol$X)
  
  n_all <- nrow(X_all)
  if (n_all == 0L) {
    stop("sobol4r_qoi_indices: design has zero rows.")
  }
  
  nrep <- as.integer(nrep)
  if (nrep < 1L) {
    stop("sobol4r_qoi_indices: nrep must be at least 1.")
  }
  
  # Evaluate the stochastic model nrep times without using replicate
  sims_list <- vector("list", nrep)
  for (i in seq_len(nrep)) {
    sims_i <- model(X_all, ...)
    if (!is.numeric(sims_i)) {
      stop("sobol4r_qoi_indices: model must return a numeric vector.")
    }
    if (length(sims_i) != n_all) {
      stop("sobol4r_qoi_indices: model output length must equal nrow(gensol$X).")
    }
    sims_list[[i]] <- sims_i
  }
  
  sims <- do.call(cbind, sims_list)
  # sims is n_all x nrep
  
  # Apply QoI function row wise
  qoi_vals <- apply(sims, 1L, qoi_fun)
  
  sensitivity::tell(gensol, qoi_vals)
}
