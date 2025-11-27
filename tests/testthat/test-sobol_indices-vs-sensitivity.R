# tests/test-sobol_indices-vs-sensitivity.R
#
# Compare Sobol4R estimators (Saltelli and Jansen) with their counterparts
# in the sensitivity package on deterministic test functions.

test_that("Saltelli estimator matches sensitivity::sobol2007 on centered output", {
  skip_if_not_installed("sensitivity")
  
  set.seed(123)
  
  n <- 1000
  d <- 3
  
  design <- sobol_design(
    n     = n,
    d     = d,
    lower = rep(-pi, d),
    upper = rep(pi, d),
    quasi = TRUE
  )
  
  A <- design$A
  B <- design$B
  
  # Reference: sensitivity::sobol2007 with manual centering
  sob_obj <- sensitivity::sobol2007(
    model = NULL,
    X1    = A,
    X2    = B,
    nboot = 0
  )
  
  Y  <- ishigami_model(sob_obj$X)
  Yc <- Y - mean(Y)
  
  sob_obj <- sensitivity::tell(sob_obj, Yc)
  
  S_sens <- as.numeric(sob_obj$S$original)
  T_sens <- as.numeric(sob_obj$T$original)
  
  # Sobol4R Saltelli implementation with internal centering
  res <- sobol_indices(
    model        = ishigami_model,
    design       = design,
    replicates   = 1L,
    estimator    = "saltelli",
    keep_samples = FALSE
  )
  
  expect_equal(res$parameters, colnames(A))
  expect_equal(res$first_order, S_sens, tolerance = 0.05)
  expect_equal(res$total_order, T_sens, tolerance = 0.05)
})

test_that("Jansen estimator matches sensitivity::soboljansen", {
  skip_if_not_installed("sensitivity")
  
  set.seed(456)
  
  n <- 1000
  d <- 3
  
  design <- sobol_design(
    n     = n,
    d     = d,
    lower = rep(-pi, d),
    upper = rep(pi, d),
    quasi = TRUE
  )
  
  A <- design$A
  B <- design$B
  
  # Reference: sensitivity::soboljansen
  sob_jan <- sensitivity::soboljansen(
    model = ishigami_model,
    X1    = A,
    X2    = B,
    nboot = 0
  )
  
  S_sens <- as.numeric(sob_jan$S$original)
  T_sens <- as.numeric(sob_jan$T$original)
  
  # Sobol4R Jansen implementation
  res <- sobol_indices(
    model        = ishigami_model,
    design       = design,
    replicates   = 1L,
    estimator    = "jansen",
    keep_samples = FALSE
  )
  
  expect_equal(res$parameters, colnames(A))
  expect_equal(res$first_order, S_sens, tolerance = 0.05)
  expect_equal(res$total_order, T_sens, tolerance = 0.05)
})
