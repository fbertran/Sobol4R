# Sobol indices for stochastic simulators with multiple estimators
# - Saltelli estimator with internal centering
# - Jansen-type estimator (robust, variance-of-differences)

#' Sobol Indices for Stochastic Simulators
#'
#' Estimate first-order and total-order Sobol indices using Monte Carlo
#' estimators that support noisy outputs via independent replicates.
#'
#' Two families of estimators are available:
#' \itemize{
#'   \item \code{"saltelli"}: Saltelli-type estimator with internal centering
#'         of the model outputs before variance and index computation.
#'   \item \code{"jansen"}: Jansen-type estimator based on variances of
#'         output differences, which is numerically stable in many settings.
#' }
#'
#' @param model Function receiving a numeric matrix and returning a numeric
#'   vector of responses. The function may include internal randomness.
#' @param design Output of [sobol_design()].
#' @param replicates Integer, number of repeated evaluations to average out the
#'   model noise. Defaults to one replicate (deterministic behaviour).
#' @param estimator Character string, either \code{"saltelli"} or
#'   \code{"jansen"}. Defaults to \code{"jansen"}.
#' @param keep_samples When \code{TRUE}, store all simulated values.
#' @param ... Further arguments passed to \code{model}.
#'
#' @return An object of class \code{sobol_result} containing the indices,
#'   intermediate estimates, and the Monte Carlo variance.
#' @export
#'
#' @examples
#' design <- sobol_design(n = 128, d = 3, quasi = TRUE)
#' model <- function(x) ishigami_model(x)
#' result <- sobol_indices(model, design, replicates = 4)
#' result$data
sobol_indices <- function(model, design, replicates = 1L,
                          estimator = c("jansen", "saltelli"),
                          keep_samples = FALSE,
                          ...) {
  estimator  <- match.arg(estimator)
  replicates <- as.integer(replicates)
  stopifnot(replicates >= 1L)
  
  design <- validate_design(design)
  A <- design$A
  B <- design$B
  
  # Evaluate model on A and B
  ya <- evaluate_model(model, A, replicates, ...)
  yb <- evaluate_model(model, B, replicates, ...)
  
  # Build mixed designs A with column j replaced by B[, j]
  mixes <- lapply(seq_len(ncol(A)), function(j) replace_column(A, B, j))
  y_mix <- lapply(
    mixes,
    evaluate_model,
    model      = model,
    replicates = replicates,
    ...
  )
  
  # Internal centering of outputs (recommended for Saltelli-type estimators)
  global_mean <- mean(c(ya$mean, yb$mean))
  ya_c        <- ya$mean - global_mean
  yb_c        <- yb$mean - global_mean
  y_mix_c     <- lapply(y_mix, function(obj) obj$mean - global_mean)
  
  # Variance of the centered output
  va <- stats::var(c(ya_c, yb_c))
  if (!is.finite(va) || va <= 0) {
    stop("Variance estimate is not positive. Increase sample size.")
  }
  
  p <- ncol(A)
  first <- numeric(p)
  total <- numeric(p)
  
  if (identical(estimator, "saltelli")) {
    # Saltelli-type estimator with centered outputs
    for (j in seq_len(p)) {
      yj <- y_mix_c[[j]]
      first[j] <- mean(yb_c * (yj - ya_c)) / va
      total[j] <- 0.5 * mean((ya_c - yj)^2) / va
    }
  } else if (identical(estimator, "jansen")) {
    # Jansen-type estimator (robust variance-of-differences)
    for (j in seq_len(p)) {
      yj <- y_mix_c[[j]]
      first[j] <- 1 - 0.5 * mean((yb_c - yj)^2) / va
      total[j] <- 0.5 * mean((ya_c - yj)^2) / va
    }
  }
  
  output <- list(
    call           = match.call(),
    parameters     = colnames(A),
    first_order    = first,
    total_order    = total,
    variance       = va,
    estimator      = estimator,
    replicates     = replicates,
    mean_A         = mean(ya$mean),
    noise_variance = mean(c(ya$variance, yb$variance)),
    data           = data.frame(
      parameter   = colnames(A),
      first_order = first,
      total_order = total
    )
  )
  
  if (isTRUE(keep_samples)) {
    output$samples <- list(A = ya, B = yb, mixes = y_mix)
  }
  
  class(output) <- "sobol_result"
  output
}

validate_design <- function(design) {
  if (!is.list(design) || !all(c("A", "B") %in% names(design))) {
    stop("design must be the output of sobol_design().")
  }
  A <- design$A
  B <- design$B
  if (!is.matrix(A) || !is.matrix(B)) {
    stop("A and B must be matrices.")
  }
  if (!all(dim(A) == dim(B))) {
    stop("A and B must have the same dimensions.")
  }
  design
}

evaluate_model <- function(model, samples, replicates, ...) {
  n <- nrow(samples)
  values <- matrix(NA_real_, nrow = n, ncol = replicates)
  for (i in seq_len(replicates)) {
    response <- model(samples, ...)
    if (length(response) != n) {
      stop("model must return a numeric vector with length equal to nrow(samples).")
    }
    values[, i] <- as.numeric(response)
  }
  list(
    mean     = rowMeans(values),
    variance = if (replicates > 1L) apply(values, 1L, stats::var) else rep(0, n),
    values   = values
  )
}

replace_column <- function(A, B, j) {
  mix <- A
  mix[, j] <- B[, j]
  mix
}