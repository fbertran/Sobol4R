#' Design generation for Sobol indices
#'
#' Simple helper that wraps \code{sensitivity::sobol} with \code{model = NULL}
#' to create the extended design matrix used to evaluate the model.
#'
#' @param X1 First sample (matrix or data.frame).
#' @param X2 Second sample (matrix or data.frame).
#' @param order Maximum interaction order (1 or 2).
#' @param nboot Number of bootstrap replicates for confidence intervals.
#' @param type Type of Monte Carlo Estimation of Sobol' Indices to be used.
#'   Supported estimators mirror the \pkg{sensitivity} helpers: \code{sobol},
#'   \code{sobol2007}, \code{soboljansen}, \code{sobolEff}, and
#'   \code{sobolmartinez}. Defaults to \code{"soboljansen"}, which is the
#'   safest general-purpose choice for both deterministic and stochastic
#'   simulators.
#' @param ... Additional arguments passed to \code{sensitivity::sobol}.
#'
#' @return An object of class \code{"sobol"} whose \code{$X} field contains
#'   the design matrix. You should evaluate your model on \code{$X} and
#'   then call \code{sensitivity::tell()}.
#' @export
sobol4r_design <- function(X1, X2, order = 2, nboot = 0, type = c("soboljansen", "sobol", "sobol2007", "sobolEff", "sobolmartinez"), ...) {
  type = match.arg(type, c("soboljansen", "sobol", "sobol2007", "sobolEff", "sobolmartinez"))
  switch(type,
         sobol = sensitivity::sobol(
           model = NULL,
           X1 = X1,
           X2 = X2,
           order = order,
           nboot = nboot,
           ...
         ),
         sobol2007 = sensitivity::sobol2007(
           model = NULL,
           X1 = X1,
           X2 = X2,
           nboot = nboot,
           ...
         ),
         soboljansen = sensitivity::soboljansen(
           model = NULL,
           X1 = X1,
           X2 = X2,
           nboot = nboot,
           ...
         ),
         sobolEff = sensitivity::sobolEff(
           model = NULL,
           X1 = X1,
           X2 = X2,
           nboot = nboot,
           ...
         ),
         sobolmartinez = sensitivity::sobolmartinez(
           model = NULL,
           X1 = X1,
           X2 = X2,
           nboot = nboot,
           ...
         )
  )
}


#' Create Sobol Sampling Designs
#'
#' Generate the two-sample matrices (A and B) that are required to apply
#' Monte Carlo Sobol estimators. The helper can rely on pseudo random numbers
#' or on a light-weight Halton low-discrepancy sequence to increase coverage.
#'
#' @param n Integer, number of rows per design matrix.
#' @param d Integer, number of model parameters.
#' @param lower Numeric vector of length d containing lower bounds.
#' @param upper Numeric vector of length d containing upper bounds.
#' @param quasi Logical, when \code{TRUE} a Halton sequence is used.
#' @param seed Optional integer used to initialise the RNG state.
#' @return A list with matrices \code{A} and \code{B} plus the column names.
#' @export
#' @examples
#' design <- sobol_design(n = 64, d = 3, quasi = TRUE)
#' str(design)
sobol_design <- function(n, d, lower = rep(0, d), upper = rep(1, d),
                         quasi = FALSE, seed = NULL) {
  stopifnot(n > 0, d > 0)
  lower <- rep(lower, length.out = d)
  upper <- rep(upper, length.out = d)
  if (!is.null(seed)) {
    set.seed(seed)
  }
  generator <- if (isTRUE(quasi)) halton_sequence else random_uniform
  if (identical(generator, halton_sequence)) {
    A <- generator(n, d, start = 0L)
    B <- generator(n, d, start = n)
  } else {
    A <- generator(n, d)
    B <- generator(n, d)
  }
  span <- upper - lower
  A <- sweep(A, 2, span, `*`)
  A <- sweep(A, 2, lower, `+`)
  B <- sweep(B, 2, span, `*`)
  B <- sweep(B, 2, lower, `+`)
  colnames(A) <- colnames(B) <- paste0("X", seq_len(d))
  list(A = A, B = B, lower = lower, upper = upper)
}

random_uniform <- function(n, d) {
  matrix(stats::runif(n * d), nrow = n, ncol = d)
}

halton_sequence <- function(n, d, start = 0L) {
  primes <- generate_primes(d)
  seqs <- lapply(primes, function(base) halton_dimension(n, base, start))
  do.call(cbind, seqs)
}

halton_dimension <- function(n, base, start = 0L) {
  result <- numeric(n)
  for (i in seq_len(n)) {
    f <- 1
    r <- i + start
    value <- 0
    while (r > 0) {
      f <- f / base
      value <- value + f * (r %% base)
      r <- floor(r / base)
    }
    result[i] <- value
  }
  result
}

generate_primes <- function(k) {
  primes <- c()
  candidate <- 2
  while (length(primes) < k) {
    if (all(candidate %% primes != 0)) {
      primes <- c(primes, candidate)
    }
    candidate <- candidate + 1
  }
  primes
}
