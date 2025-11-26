#' Estimate Failure Probability from Simulator Outputs
#'
#' Convenient helper to compute the reliability-related probabilities described
#' in Lebrun et al. (2021). The failure domain is controlled by a threshold and
#' an inequality direction.
#'
#' @param response Numeric vector of simulator evaluations.
#' @param threshold Numeric scalar defining the failure boundary.
#' @param less Logical, failure is defined as \code{response <= threshold} when
#'   \code{TRUE} and \code{response >= threshold} otherwise.
#' @param weights Optional numeric vector of non-negative weights. The vector is
#'   normalised internally when supplied.
#' @return A list containing the estimated probability and its variance.
#' @export
#' @examples
#' y <- rnorm(1000)
#' estimate_failure_probability(y, threshold = -1)
estimate_failure_probability <- function(response, threshold, less = TRUE,
                                         weights = NULL) {
  response <- as.numeric(response)
  if (length(threshold) != 1L || !is.finite(threshold)) {
    stop("threshold must be a single finite numeric value.")
  }
  valid <- is.finite(response)
  if (!all(valid)) {
    response <- response[valid]
    if (!is.null(weights)) {
      weights <- weights[valid]
    }
  }
  if (length(response) == 0L) {
    stop("response must contain at least one finite value.")
  }
  if (is.null(weights)) {
    weights <- rep(1 / length(response), length(response))
  } else {
    if (length(weights) != length(response)) {
      stop("weights must match the length of response.")
    }
    if (any(!is.finite(weights)) || any(weights < 0)) {
      stop("weights must be non-negative and finite.")
    }
    weights <- weights / sum(weights)
  }
  indicator <- if (isTRUE(less)) response <= threshold else response >= threshold
  probability <- sum(weights * indicator)
  variance <- probability * (1 - probability) * sum(weights^2)
  list(probability = probability, variance = variance)
}