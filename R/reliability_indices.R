#' Reliability-Oriented Sobol Indices
#'
#' Transform stored simulator samples into Sobol indices for the binary failure
#' indicator described by Lebrun et al. (2021). The function reuses the
#' Saltelli-type estimator from [sobol_indices()] and therefore requires a
#' previous call with `keep_samples = TRUE`.
#'
#' @param result Output of [sobol_indices()] computed with
#'   `keep_samples = TRUE`.
#' @param threshold Numeric scalar defining the failure boundary.
#' @param less Logical, when `TRUE` failures correspond to
#'   `response <= threshold`; otherwise, failures correspond to
#'   `response >= threshold`.
#' @return A `sobol_result` instance storing the Sobol indices of the failure
#'   indicator along with the estimated failure probability and its variance.
#' @examples
#' design <- sobol_design(n = 128, d = 3, lower = rep(-pi, 3), upper = rep(pi, 3))
#' stochastic <- sobol_indices(ishigami_model, design, replicates = 3,
#'                             keep_samples = TRUE)
#' failure <- sobol_reliability(stochastic, threshold = -1)
#' Sobol4R::autoplot(failure, show_uncertainty = TRUE)
#' @export
sobol_reliability <- function(result, threshold, less = TRUE) {
  stopifnot(inherits(result, "sobol_result"))
  samples <- result$samples
  if (is.null(samples)) {
    stop("sobol_reliability() requires 'keep_samples = TRUE' in sobol_indices().")
  }
  if (length(threshold) != 1L || !is.finite(threshold)) {
    stop("threshold must be a single finite numeric value.")
  }
  parameters <- result$parameters
  if (is.null(parameters)) {
    parameters <- names(samples$mixes)
  }
  parameters <- result$parameters
  if (is.null(parameters)) {
    parameters <- paste0("X", seq_along(samples$mixes))
  }
  indicator_A <- failure_indicator(samples$A, threshold, less)
  indicator_B <- failure_indicator(samples$B, threshold, less)
  indicator_mixes <- lapply(samples$mixes, failure_indicator,
                            threshold = threshold, less = less)
  va <- stats::var(c(indicator_A, indicator_B))
  if (!is.finite(va) || va <= 0) {
    stop("Failure indicator variance is zero; adjust the threshold or sample size.")
  }
  first <- vapply(indicator_mixes, function(mix) {
    mean(indicator_B * (mix - indicator_A)) / va
  }, numeric(1L))
  total <- vapply(indicator_mixes, function(mix) {
    0.5 * mean((indicator_A - mix)^2) / va
  }, numeric(1L))
  weights_pb <- rep(1 / (length(indicator_A) + length(indicator_B)),
                    length(indicator_A) + length(indicator_B))
  failure_stats <- estimate_failure_probability(
    response = c(indicator_A, indicator_B),
    threshold = 0.5,
    less = FALSE,
    weights = weights_pb
  )
  output <- list(
    call = match.call(),
    parameters = parameters,
    first_order = as.numeric(first),
    total_order = as.numeric(total),
    variance = va,
    estimator = "saltelli",
    replicates = 1L,
    mean_A = mean(indicator_A),
    noise_variance = 0,
    failure_probability = failure_stats$probability,
    failure_variance = failure_stats$variance,
    threshold = threshold,
    less = less,
    data = data.frame(parameter = parameters,
                      first_order = as.numeric(first),
                      total_order = as.numeric(total))
  )
  output$samples <- list(
    A = list(values = matrix(indicator_A, ncol = 1L)),
    B = list(values = matrix(indicator_B, ncol = 1L)),
    mixes = lapply(indicator_mixes, function(mix) {
      list(values = matrix(mix, ncol = 1L))
    })
  )
  class(output) <- "sobol_result"
  output
}

failure_indicator <- function(block, threshold, less) {
  values <- block$values
  if (is.null(values)) {
    stop("Stored samples are missing raw 'values'. Re-run sobol_indices() with keep_samples = TRUE.")
  }
  values <- as.matrix(values)
  if (!is.numeric(values)) {
    stop("Stored samples must be numeric.")
  }
  averages <- rowMeans(values, na.rm = TRUE)
  if (any(is.na(averages))) {
    stop("Failure indicator encountered rows with only missing values.")
  }
  indicator <- if (isTRUE(less)) averages <= threshold else averages >= threshold
  as.numeric(indicator)
}