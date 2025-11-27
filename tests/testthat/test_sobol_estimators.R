skip_if_not_installed("sensitivity")

rng_reset <- function() {
  set.seed(42)
}

test_that("sobol4r_design accepts all supported estimators", {
  rng_reset()
  X1 <- matrix(runif(96), ncol = 3)
  X2 <- matrix(runif(96), ncol = 3)
  types <- c("sobol", "sobol2007", "soboljansen", "sobolEff", "sobolmartinez")
  designs <- lapply(types, function(t) sobol4r_design(X1 = X1, X2 = X2, type = t, nboot = 0))
  mapply(function(obj, cls) expect_s3_class(obj, cls), designs, types)
})

test_that("sobol4r_run returns filled objects for each estimator", {
  rng_reset()
  n <- 64
  d <- 3
  X1 <- matrix(runif(n * d), ncol = d)
  X2 <- matrix(runif(n * d), ncol = d)
  types <- c("sobol", "sobol2007", "soboljansen", "sobolEff", "sobolmartinez")
  model <- function(X) ishigami_model(as.matrix(X))
  
  for (typ in types) {
    res <- sobol4r_run(model, X1 = X1, X2 = X2, type = typ, nboot = 0)
    expect_s3_class(res, typ)
    if (!inherits(res, "sobol") && !inherits(res, "sobolEff")) {
      expect_true(!is.null(res$T))
    }
    expect_true(!is.null(res$S))
  }
})

test_that("autoplot supports all sensitivity Sobol classes", {
  skip_if_not_installed("ggplot2")
  rng_reset()
  X1 <- matrix(runif(48), ncol = 3)
  X2 <- matrix(runif(48), ncol = 3)
  model <- function(X) ishigami_model(as.matrix(X))
  types <- c("sobol2007", "soboljansen", "sobolmartinez")
  for (t in types) {
    res <- sobol4r_run(model, X1 = X1, X2 = X2, type = t, nboot = 100)
    plt <- autoplot(res)
    expect_s3_class(plt, "ggplot")
  }
})
