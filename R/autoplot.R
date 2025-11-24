#' Summarise Sobol Indices
#'
#' Compute compact summaries of the Sobol indices and their Monte Carlo
#' variability. The function is intended to feed diagnostic plots.
#'
#' @param result A \code{sobol_result} object.
#' @param probs Numeric vector of probabilities used to report quantiles of the
#'   empirical bootstrap distribution.
#' @param bootstrap Integer, number of bootstrap resamples used to quantify
#'   the estimator uncertainty.
#' @return A data frame (class \code{sobol_summary}) with the requested
#'   statistics. Quantile columns are added when \code{probs} is not empty.
#' @export
#' @examples
#' design <- sobol_design(n = 64, d = 3)
#' model <- function(x) ishigami_model(x)
#' sob <- sobol_indices(model, design, keep_samples = TRUE)
#' summarise_sobol(sob, probs = c(0.1, 0.9))
summarise_sobol <- function(result, probs = c(0.1, 0.5, 0.9), bootstrap = 200L) {
  stopifnot(inherits(result, "sobol_result"))
  probs <- sort(unique(probs))
  if (length(probs) && any(probs < 0 | probs > 1)) {
    stop("probs must lie within [0, 1].")
  }
  bootstrap <- as.integer(bootstrap)
  if (length(bootstrap) != 1L || is.na(bootstrap) || bootstrap < 0L) {
    stop("bootstrap must be a non-negative integer.")
  }
  stats <- bootstrap_indices(result, bootstrap)
  summary <- data.frame(
    parameter = result$parameters,
    first_order = result$first_order,
    total_order = result$total_order,
    variance = result$variance,
    noise_variance = result$noise_variance,
    row.names = NULL, check.names = FALSE
  )
  if (length(probs)) {
    quantiles_first <- t(vapply(seq_len(ncol(stats$first)), function(j) {
      stats::quantile(stats$first[, j], probs = probs, na.rm = TRUE, names = FALSE)
    }, numeric(length(probs))))
    quantiles_total <- t(vapply(seq_len(ncol(stats$total)), function(j) {
      stats::quantile(stats$total[, j], probs = probs, na.rm = TRUE, names = FALSE)
    }, numeric(length(probs))))
    prob_labels <- format_probabilities(probs)
    for (k in seq_along(prob_labels)) {
      summary[[paste0("first_q", prob_labels[k])]] <- quantiles_first[, k]
      summary[[paste0("total_q", prob_labels[k])]] <- quantiles_total[, k]
    }
  }
  class(summary) <- c("sobol_summary", class(summary))
  attr(summary, "probs") <- probs
  summary
}

#' @export
bootstrap_indices <- function(result, bootstrap) {
  p <- length(result$parameters)
  if (is.null(result$samples) || bootstrap < 1L) {
    return(list(first = matrix(result$first_order, nrow = 1, ncol = p, byrow = TRUE),
                total = matrix(result$total_order, nrow = 1, ncol = p, byrow = TRUE)))
  }
  A <- result$samples$A$values
  mixes <- lapply(result$samples$mixes, function(x) x$values)
  n <- nrow(A)
  first_boot <- matrix(NA_real_, nrow = bootstrap, ncol = length(mixes))
  total_boot <- matrix(NA_real_, nrow = bootstrap, ncol = length(mixes))
  for (b in seq_len(bootstrap)) {
    idx <- sample.int(n, replace = TRUE)
    ya <- rowMeans(result$samples$A$values[idx, , drop = FALSE])
    yb <- rowMeans(result$samples$B$values[idx, , drop = FALSE])
    va <- stats::var(c(ya, yb))
    for (j in seq_along(mixes)) {
      y_mix <- rowMeans(mixes[[j]][idx, , drop = FALSE])
      first_boot[b, j] <- mean(yb * (y_mix - ya)) / va
      total_boot[b, j] <- 0.5 * mean((ya - y_mix)^2) / va
    }
  }
  list(first = first_boot, total = total_boot)
}

format_probabilities <- function(probs) {
  labels <- gsub("\\.", "_", format(probs, trim = TRUE, scientific = FALSE))
  labels <- sub("_0+$", "", labels)
  labels <- sub("^0_", "0", labels)
  labels
}

#' Generic autoplot implementation
#'
#' Provide a ggplot visualisation when ggplot2 is available, otherwise fallback
#' to a lightweight base R bar chart. Supports the custom `sobol_result` class
#' used in this package, compact `sobol_summary` data frames, and
#' `sensitivity::sobol` objects.
#'
#' @param ' @param object A \code{sobol_result}, \code{sobol_summary}, or
#'   \code{sensitivity::sobol} instance.
#' @param show_uncertainty Logical, when \code{TRUE} bootstrap quantiles are
#'   computed (if available) and displayed as error bars.
#' @param probs Numeric vector of probabilities used for the uncertainty bars.
#' @param bootstrap Integer indicating how many bootstrap resamples to draw when
#'   \code{show_uncertainty = TRUE}.
#' @param ... Further arguments passed to the plotting backend.
#' @return A ggplot object when \code{ggplot2} is installed, otherwise the
#'   bar centres invisibly.
#' @export
autoplot <- function(object, ...) {
  UseMethod("autoplot")
}

#' @export
autoplot.sobol_result <- function(object, show_uncertainty = FALSE,
                                  probs = c(0.1, 0.9), bootstrap = 200L, ...) {
  if (isTRUE(show_uncertainty)) {
    summary <- summarise_sobol(object, probs = probs, bootstrap = bootstrap)
    return(autoplot(summary, ...))
  }
  data <- object$data
  if (is.null(data)) {
    parameters <- object$parameters
    if (is.null(parameters)) {
      parameters <- paste0("X", seq_along(object$first_order))
    }
    data <- data.frame(
      parameter = parameters,
      first_order = object$first_order,
      total_order = object$total_order,
      row.names = NULL, check.names = FALSE
    )
  }
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    df <- data.frame(
      parameter = rep(data$parameter, 2),
      index_type = rep(c("First", "Total"), each = nrow(data)),
      value = c(data$first_order, data$total_order)
    )
    p <- ggplot2::ggplot(df, ggplot2::aes(x = parameter, y = value,
                                          fill = index_type)) +
      ggplot2::geom_col(position = "dodge") +
      ggplot2::geom_hline(yintercept = 0, colour = "grey50") +
      ggplot2::labs(y = "Sobol index", fill = "Type",
                    title = "Sobol sensitivity summary") +
      ggplot2::theme_minimal()
    return(p)
  }
  graphics::barplot(rbind(data$first_order, data$total_order),
                    beside = TRUE, names.arg = data$parameter,
                    legend.text = c("First", "Total"),
                    args.legend = list(x = "topright"), ...)
  invisible(NULL)
}

#' @export
autoplot.sobol <- function(object, ...) {
  stopifnot(inherits(object, "sobol"))
  S <- object$S
  T <- object$T
  if (is.null(S)) {
    stop("sobol object is missing the 'S'. Run sensitivity::tell().")
  }
  extract_indices <- function(df) {
    if (!is.data.frame(df)) {
      stop("Sobol index tables must be data frames.")
    }
    parameters <- rownames(df)
    if (is.null(parameters) || any(parameters == "")) {
      parameters <- seq_len(nrow(df))
    }
    values <- if ("original" %in% names(df)) df$original else df[[1L]]
    interval <- NULL
    if (all(c("min. c.i.", "max. c.i.") %in% names(df))) {
      interval <- cbind(df[["min. c.i."]], df[["max. c.i."]])
    }
    list(parameters = parameters, values = values, interval = interval)
  }
  first_info <- extract_indices(S)
  total_info <- extract_indices(T)
  df <- data.frame(
    parameter = rep(first_info$parameters, 2),
    index_type = rep(c("First", "Total"), each = length(first_info$parameters)),
    value = c(first_info$values, total_info$values)
  )
  lower <- upper <- rep(NA_real_, nrow(df))
  if (!is.null(first_info$interval)) {
    lower[seq_len(nrow(S))] <- first_info$interval[, 1]
    upper[seq_len(nrow(S))] <- first_info$interval[, 2]
  }
  if (!is.null(total_info$interval)) {
    idx <- seq_len(nrow(T)) + nrow(S)
    lower[idx] <- total_info$interval[, 1]
    upper[idx] <- total_info$interval[, 2]
  }
  df$lower <- lower
  df$upper <- upper
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = parameter, y = value,
                                          fill = index_type)) +
      ggplot2::geom_col(position = "dodge") +
      ggplot2::geom_hline(yintercept = 0, colour = "grey50") +
      ggplot2::labs(y = "Sobol index", fill = "Type",
                    title = "Sobol sensitivity summary") +
      ggplot2::theme_minimal()
    if (any(!is.na(df$lower))) {
      p <- p + ggplot2::geom_errorbar(
        data = df,
        ggplot2::aes(x = parameter, ymin = lower, ymax = upper,
                     group = index_type),
        position = ggplot2::position_dodge(width = 0.9), width = 0.2,
        size = 0.4, inherit.aes = FALSE)
    }
    return(p)
  }
  heights <- rbind(first_info$values, total_info$values)
  bars <- graphics::barplot(heights, beside = TRUE,
                            names.arg = first_info$parameters,
                            legend.text = c("First", "Total"),
                            args.legend = list(x = "topright"), ...)
  if (!is.null(first_info$interval)) {
    graphics::arrows(x0 = bars[1, ], x1 = bars[1, ],
                     y0 = first_info$interval[, 1], y1 = first_info$interval[, 2],
                     angle = 90, code = 3, length = 0.05)
  }
  if (!is.null(total_info$interval)) {
    graphics::arrows(x0 = bars[2, ], x1 = bars[2, ],
                     y0 = total_info$interval[, 1], y1 = total_info$interval[, 2],
                     angle = 90, code = 3, length = 0.05)
  }
  invisible(NULL)
}

#' @export
autoplot.sobol2007 <- function(object, ...) {
  stopifnot(inherits(object, "sobol2007"))
  S <- object$S
  T <- object$T
  if (is.null(S) || is.null(T)) {
    stop("sobol object is missing the 'S' or 'T' components. Run sensitivity::tell().")
  }
  extract_indices <- function(df) {
    if (!is.data.frame(df)) {
      stop("Sobol index tables must be data frames.")
    }
    parameters <- rownames(df)
    if (is.null(parameters) || any(parameters == "")) {
      parameters <- seq_len(nrow(df))
    }
    values <- if ("original" %in% names(df)) df$original else df[[1L]]
    interval <- NULL
    if (all(c("min. c.i.", "max. c.i.") %in% names(df))) {
      interval <- cbind(df[["min. c.i."]], df[["max. c.i."]])
    }
    list(parameters = parameters, values = values, interval = interval)
  }
  first_info <- extract_indices(S)
  total_info <- extract_indices(T)
  df <- data.frame(
    parameter = rep(first_info$parameters, 2),
    index_type = rep(c("First", "Total"), each = length(first_info$parameters)),
    value = c(first_info$values, total_info$values)
  )
  lower <- upper <- rep(NA_real_, nrow(df))
  if (!is.null(first_info$interval)) {
    lower[seq_len(nrow(S))] <- first_info$interval[, 1]
    upper[seq_len(nrow(S))] <- first_info$interval[, 2]
  }
  if (!is.null(total_info$interval)) {
    idx <- seq_len(nrow(T)) + nrow(S)
    lower[idx] <- total_info$interval[, 1]
    upper[idx] <- total_info$interval[, 2]
  }
  df$lower <- lower
  df$upper <- upper
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = parameter, y = value,
                                          fill = index_type)) +
      ggplot2::geom_col(position = "dodge") +
      ggplot2::geom_hline(yintercept = 0, colour = "grey50") +
      ggplot2::labs(y = "Sobol index", fill = "Type",
                    title = "Sobol sensitivity summary") +
      ggplot2::theme_minimal()
    if (any(!is.na(df$lower))) {
      p <- p + ggplot2::geom_errorbar(
        data = df,
        ggplot2::aes(x = parameter, ymin = lower, ymax = upper,
                     group = index_type),
        position = ggplot2::position_dodge(width = 0.9), width = 0.2,
        size = 0.4, inherit.aes = FALSE)
    }
    return(p)
  }
  heights <- rbind(first_info$values, total_info$values)
  bars <- graphics::barplot(heights, beside = TRUE,
                            names.arg = first_info$parameters,
                            legend.text = c("First", "Total"),
                            args.legend = list(x = "topright"), ...)
  if (!is.null(first_info$interval)) {
    graphics::arrows(x0 = bars[1, ], x1 = bars[1, ],
                     y0 = first_info$interval[, 1], y1 = first_info$interval[, 2],
                     angle = 90, code = 3, length = 0.05)
  }
  if (!is.null(total_info$interval)) {
    graphics::arrows(x0 = bars[2, ], x1 = bars[2, ],
                     y0 = total_info$interval[, 1], y1 = total_info$interval[, 2],
                     angle = 90, code = 3, length = 0.05)
  }
  invisible(NULL)
}

#' @export
autoplot.sobol_summary <- function(object, ...) {
  stopifnot(inherits(object, "sobol_summary"))
  quant_cols_first <- grep("^first_q", names(object), value = TRUE)
  quant_cols_total <- grep("^total_q", names(object), value = TRUE)
  interval_first <- interval_total <- NULL
  if (length(quant_cols_first) >= 2L) {
    interval_first <- cbind(object[[quant_cols_first[1L]]],
                            object[[quant_cols_first[length(quant_cols_first)]]])
  }
  if (length(quant_cols_total) >= 2L) {
    interval_total <- cbind(object[[quant_cols_total[1L]]],
                            object[[quant_cols_total[length(quant_cols_total)]]])
  }
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    df <- data.frame(
      parameter = rep(object$parameter, 2),
      index_type = rep(c("First", "Total"), each = nrow(object)),
      value = c(object$first_order, object$total_order)
    )
    lower <- upper <- rep(NA_real_, nrow(df))
    if (!is.null(interval_first)) {
      lower[seq_len(nrow(object))] <- interval_first[, 1]
      upper[seq_len(nrow(object))] <- interval_first[, 2]
    }
    if (!is.null(interval_total)) {
      idx <- seq_len(nrow(object)) + nrow(object)
      lower[idx] <- interval_total[, 1]
      upper[idx] <- interval_total[, 2]
    }
    df$lower <- lower
    df$upper <- upper
    p <- ggplot2::ggplot(df, ggplot2::aes(x = parameter, y = value,
                                          fill = index_type)) +
      ggplot2::geom_col(position = "dodge") +
      ggplot2::geom_hline(yintercept = 0, colour = "grey50") +
      ggplot2::labs(y = "Sobol index", fill = "Type",
                    title = "Sobol sensitivity summary") +
      ggplot2::theme_minimal()
    if (any(!is.na(df$lower))) {
      p <- p + ggplot2::geom_errorbar(
        data = df,
        ggplot2::aes(x = parameter, ymin = lower, ymax = upper,
                     group = index_type),
        position = ggplot2::position_dodge(width = 0.9), width = 0.2,
        size = 0.4, inherit.aes = FALSE)
    }
    return(p)
  }
  heights <- rbind(object$first_order, object$total_order)
  bars <- graphics::barplot(heights, beside = TRUE, names.arg = object$parameter,
                            legend.text = c("First", "Total"),
                            args.legend = list(x = "topright"), ...)
  if (!is.null(interval_first)) {
    graphics::arrows(x0 = bars[1, ], x1 = bars[1, ],
                     y0 = interval_first[, 1], y1 = interval_first[, 2],
                     angle = 90, code = 3, length = 0.05)
  }
  if (!is.null(interval_total)) {
    graphics::arrows(x0 = bars[2, ], x1 = bars[2, ],
                     y0 = interval_total[, 1], y1 = interval_total[, 2],
                     angle = 90, code = 3, length = 0.05)
  }
  invisible(NULL)
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
  A <- generator(n, d)
  B <- generator(n, d)
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

halton_sequence <- function(n, d) {
  primes <- generate_primes(d)
  seqs <- lapply(primes, function(base) halton_dimension(n, base))
  do.call(cbind, seqs)
}

halton_dimension <- function(n, base) {
  result <- numeric(n)
  for (i in seq_len(n)) {
    f <- 1
    r <- i
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