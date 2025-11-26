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

#' Bootstrap Sobol indices from stored samples
#'
#' Recompute Sobol first- and total-order indices from stored sample matrices
#' using bootstrap resampling. Falls back to deterministic values when no
#' samples are available.
#'
#' @param result A \code{sobol_result} object produced by
#'   \code{sobol_indices()}.
#' @param bootstrap Integer indicating how many bootstrap replicates to draw.
#' @return A list with matrices \code{first} and \code{total} containing the
#'   bootstrap replications.
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

#' Autoplot implementations
#'
#' Provide a ggplot visualisation when ggplot2 is available, otherwise fallback
#' to a lightweight base R bar chart. Supports the custom `sobol_result` class
#' used in this package, compact `sobol_summary` data frames, and
#' `sensitivity::sobol` objects.
#'
#' @name Autoplot implementations
#' @rdname autoplot
#' @param object A \code{sobol_result}, \code{sobol_summary}, or
#'   \code{sensitivity::sobol} instance.
#' @param show_uncertainty Logical, when \code{TRUE} bootstrap quantiles are
#'   computed (if available) and displayed as error bars.
#' @param probs Numeric vector of probabilities used for the uncertainty bars.
#' @param bootstrap Integer indicating how many bootstrap resamples to draw when
#'   \code{show_uncertainty = TRUE}.
#' @param separate_panels Should the indices be plotted on separate 
#'   panels according to their order? 
#'   If `separate_panels = TRUE`, the first order indices are separated from 
#'   the higher orders ones.
#' @param ncol If `separate_panels = TRUE`, the number of columns for the 
#'   facet wrapping of the plot.
#' @param ... Further arguments passed to the plotting backend.
#' @return A ggplot object when \code{ggplot2} is installed, otherwise the
#'   bar centres invisibly.
NULL

#' Generic autoplot for `Sobol4R` package
#' 
#' @rdname autoplot
#' @export
autoplot <- function(object, ...) {
  UseMethod("autoplot")
}

#' Autoplot for `sobol_result` class
#' 
#' @rdname autoplot
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
  if (rlang::is_installed("ggplot2")) {
    df <- data.frame(
      parameter = rep(data$parameter, 2),
      index_type = rep(c("First", "Total"), each = nrow(data)),
      value = c(data$first_order, data$total_order)
    )
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$parameter, y = .data$value,
                                          fill = .data$index_type)) +
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

#' Autoplot for `sobol` class
#' 
#' @rdname autoplot
#' @export
autoplot.sobol <- function(object, separate_panels = TRUE, ncol= 2, ...) {
  stopifnot(inherits(object, "sobol"))
  
  S  <- object$S        # first order
  
  if (is.null(S)) stop("sobol object missing S. Run tell().")
  # object$T is always NULL for sensitivity::sobol()
  
  # Accept both X1:X2 and X1*X2 formats
  # Helper: parse "X1:X2" -> degree 2
  parse_degree <- function(name) {
    name2 <- gsub("\\*", ":", name)  # convert X1*X2 → X1:X2 internally
    vars <- unlist(strsplit(name2, ":"))
    vars <- vars[vars != ""]
    length(vars)
  }
  
  ind_deg1 <- vapply(rownames(S), parse_degree, numeric(1))==1
  Sh <- object$Shigh    # higher order (may contain 2nd, 3rd, ...)
  Sh <- rbind(Sh, S[!ind_deg1,])
  
  S <- S[ind_deg1,]
  
  # Helper: convert a Sobol index data.frame into a clean long format
  extract_df <- function(tab, label) {
    if (is.null(tab)) return(NULL)
    if (!is.data.frame(tab)) stop("Sobol index tables must be data frames.")
    
    par <- rownames(tab)
    if (is.null(par)) par <- seq_len(nrow(tab))
    
    # Replace ":" with "*" in rownames for display
    par_display <- gsub(":", "*", par)
    
    val <- if ("original" %in% colnames(tab)) tab$original else tab[[1]]
    degree <- vapply(par_display, parse_degree, numeric(1))
    
    ci_low  <- if ("min. c.i." %in% names(tab)) tab[["min. c.i."]] else NA_real_
    ci_high <- if ("max. c.i." %in% names(tab)) tab[["max. c.i."]] else NA_real_
    
    data.frame(
      parameter  = par_display,
      degree     = degree,
      index_type = label,
      value      = val,
      lower      = ci_low,
      upper      = ci_high,
      stringsAsFactors = FALSE
    )
  }
  
  # Build dataset — note: no "Total Order" group for sensitivity::sobol()
  if(nrow(Sh)>0){
  df <- do.call(
    rbind,
    list(
      extract_df(S,  "First order"),
      extract_df(Sh, "Higher order")
    )
  )} else{
    df <- do.call(
      rbind,
      list(
        extract_df(S,  "First order")
      )
    )
  }
  
  # Pretty degree labels
  unique_deg <- sort(unique(df$degree))
  df$degree <- factor(df$degree,
                      levels = unique_deg,
                      labels = paste("Degree", unique_deg))
  
  # ---- ggplot version ----
  if (rlang::is_installed("ggplot2")) {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$parameter, y = .data$value, fill = .data$index_type)) +
      ggplot2::geom_col(position = "dodge") +
      ggplot2::geom_hline(yintercept = 0, color = "grey50") +
      ggplot2::labs(
        title = "Sobol' Sensitivity Indices (grouped by interaction order)",
        y     = "Sobol index",
        fill  = "Index type"
      ) +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(
        angle = 45, hjust = 1)
      )
    
    # Confidence intervals
    if (any(!is.na(df$lower))) {
      p <- p +
        ggplot2::geom_errorbar(
          ggplot2::aes(ymin = .data$lower, ymax = .data$upper, group = .data$index_type),          width     = 0.2,
          linewidth      = 0.4,
          position  = ggplot2::position_dodge(width = 0.9)
        )
    }
    
    # Facet per degree (1st order, 2nd order, etc.)
    if (separate_panels) {
      p <- p + ggplot2::facet_wrap(~index_type, scales = "free_x", dir="h", ncol=ncol)
    }
    
    return(p)
  }
  
  # ---- base R fallback ----
  split_deg <- split(df, df$degree)
  for (nm in names(split_deg)) {
    message("Plotting ", nm)
    d <- split_deg[[nm]]
    graphics::barplot(d$value, names.arg = d$parameter, main = nm, ...)
  }
  
  invisible(df)
}


#' Autoplot for `sobol2007` class
#' 
#' @rdname autoplot
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
  if (rlang::is_installed("ggplot2")) {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$parameter, y = .data$value,
                                          fill = .data$index_type)) +
      ggplot2::geom_col(position = "dodge") +
      ggplot2::geom_hline(yintercept = 0, colour = "grey50") +
      ggplot2::labs(y = "Sobol index", fill = "Type",
                    title = "Sobol sensitivity summary") +
      ggplot2::theme_minimal()
    if (any(!is.na(df$lower))) {
      p <- p + ggplot2::geom_errorbar(
        data = df,
        ggplot2::aes(x = .data$parameter, ymin = .data$lower, ymax = .data$upper,
                     group = .data$index_type),
        position = ggplot2::position_dodge(width = 0.9), width = 0.2,
        linewidth = 0.4, inherit.aes = FALSE)
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

#' Autoplot for `sobol_summary` class
#' 
#' @rdname autoplot
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
  if (rlang::is_installed("ggplot2")) {
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
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$parameter, y = .data$value,
                                          fill = .data$index_type)) +
      ggplot2::geom_col(position = "dodge") +
      ggplot2::geom_hline(yintercept = 0, colour = "grey50") +
      ggplot2::labs(y = "Sobol index", fill = "Type",
                    title = "Sobol sensitivity summary") +
      ggplot2::theme_minimal()
    if (any(!is.na(df$lower))) {
      p <- p + ggplot2::geom_errorbar(
        data = df,
        ggplot2::aes(x = .data$parameter, ymin = .data$lower, ymax = .data$upper,
                     group = .data$index_type),
        position = ggplot2::position_dodge(width = 0.9), width = 0.2,
        linewidth = 0.4, inherit.aes = FALSE)
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