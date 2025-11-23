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
#' to a lightweight base R bar chart.
#'
#' @param object A \code{sobol_result} or \code{sobol_summary} instance.
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
  stopifnot(is.data.frame(data))
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
