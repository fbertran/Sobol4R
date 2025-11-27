test_that("ishigami Sobol indices agree across estimators when using the same domain", {
  set.seed(123)

  design <- sobol_design(
    n     = 512,
    d     = 3,
    lower = rep(-pi, 3),
    upper = rep(pi, 3),
    quasi = TRUE
  )

  sobol4r_result <- sobol_indices(ishigami_model, design, replicates = 1, estimator = "saltelli")

  sens_design_s7 <- sobol4r_design(
    X1 = as.data.frame(design$A),
    X2 = as.data.frame(design$B),
    order = 1,
    nboot = 0,
    type = "sobol2007"
  )
  sens_output_s7 <- sensitivity::tell(sens_design_s7, ishigami_model(as.matrix(sens_design_s7$X)))

  expect_equal(
    sobol4r_result$data$first_order,
    sens_output_s7$S$original,
    tolerance = 0.05
  )
  expect_equal(
    sobol4r_result$data$total_order,
    sens_output_s7$T$original,
    tolerance = 0.05
  )
  
  sens_design <- sobol4r_design(
    X1 = as.data.frame(design$A),
    X2 = as.data.frame(design$B),
    order = 1,
    nboot = 0,
    type = "sobol"
  )
  Y <- ishigami_model(as.matrix(sens_design$X))
  Yc <- Y - mean(Y)
  sens_output <- sensitivity::tell(sens_design, Yc)
  
})

