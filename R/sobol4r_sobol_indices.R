#' Sobol Indices for Stochastic Simulators
#'
#' Estimate first-order and total-order Sobol indices using a Saltelli-type
#' Monte Carlo estimator that supports noisy outputs via independent replicates.
#'
#' @param model Function receiving a numeric matrix and returning a numeric
#'   vector of responses. The function may include internal randomness.
#' @param design Output of [sobol_design()].
#' @param replicates Integer, number of repeated evaluations to average out the
#'   model noise. Defaults to one replicate (deterministic behaviour).
#' @param estimator Currently only \code{"saltelli"} is implemented.
#' @param keep_samples When \code{TRUE}, store all simulated values.
#' @param ... Further arguments passed to \code{model}.
#' @return An object of class \code{sobol_result} containing the indices,
#'   intermediate estimates, and the Monte Carlo variance.
#' @export
#' @examples
#' design <- sobol_design(n = 128, d = 3, quasi = TRUE)
#' model <- function(x) ishigami_model(x)
#' result <- sobol_indices(model, design, replicates = 4)
#' autoplot(result)
sobol_indices <- function(model, design, replicates = 1L,
                          estimator = c("saltelli"), keep_samples = FALSE,
                          ...) {
  estimator <- match.arg(estimator)
  replicates <- as.integer(replicates)
  stopifnot(replicates >= 1L)
  design <- validate_design(design)
  A <- design$A
  B <- design$B
  ya <- evaluate_model(model, A, replicates, ...)
  yb <- evaluate_model(model, B, replicates, ...)
  mixes <- lapply(seq_len(ncol(A)), function(j) replace_column(A, B, j))
  y_mix <- lapply(mixes, evaluate_model, model = model, replicates = replicates,
                  ...)
  va <- stats::var(c(ya$mean, yb$mean))
  if (!is.finite(va) || va <= 0) {
    stop("Variance estimate is not positive. Increase sample size.")
  }
  first <- numeric(ncol(A))
  total <- numeric(ncol(A))
  for (j in seq_along(y_mix)) {
    y_mix_mean <- y_mix[[j]]$mean
    first[j] <- mean(yb$mean * (y_mix_mean - ya$mean)) / va
    total[j] <- 0.5 * mean((ya$mean - y_mix_mean)^2) / va
  }
  output <- list(
    call = match.call(),
    parameters = colnames(A),
    first_order = first,
    total_order = total,
    variance = va,
    estimator = estimator,
    replicates = replicates,
    mean_A = mean(ya$mean),
    noise_variance = mean(c(ya$variance, yb$variance)),
    data = data.frame(parameter = colnames(A),
                      first_order = first,
                      total_order = total)
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
  list(mean = rowMeans(values),
       variance = if (replicates > 1L) apply(values, 1, stats::var) else rep(0, n),
       values = values)
}

replace_column <- function(A, B, j) {
  mix <- A
  mix[, j] <- B[, j]
  mix
}