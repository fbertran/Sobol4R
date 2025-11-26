test_that("sobol_g_cpp matches pure R implementation", {
  set.seed(123)
  X <- matrix(runif(5 * 4), ncol = 4)
  a <- c(0, 1, 4.5, 9, 99, 99, 99, 99)
  
  g_R   <- sobol_g_R(X, a = a)
  g_cpp <- sobol_g_function(X, a = a)
  
  expect_equal(g_cpp, g_R)
})

test_that("sobol_g2_cpp matches pure R implementation", {
  set.seed(123)
  X <- matrix(runif(10 * 3), ncol = 3)
  a <- c(0, 1, 4.5, 9, 99, 99, 99, 99)
  
  g2_R   <- sobol_g2_R(X, a = a)
  g2_cpp <- sobol_g2_function(X, a = a)
  
  expect_equal(g2_cpp, g2_R)
})

test_that("sobol_g2_additive_noise_cpp reproduces R version under same seed", {
  set.seed(123)
  X <- matrix(runif(20 * 2), ncol = 2)
  a <- c(0, 1, 4.5, 9, 99, 99, 99, 99)
  sd <- 0.7
  
  set.seed(42)
  y_R <- sobol_g2_additive_noise_R(X, sd = sd, a = a)
  set.seed(42)
  y_cpp <- sobol_g2_additive_noise(X, sd = sd, a = a)
  
  expect_equal(y_cpp, y_R)
})

test_that("sobol_g2_qoi_mean_cpp reproduces R version under same seed", {
  set.seed(123)
  X <- matrix(runif(15 * 2), ncol = 2)
  a <- c(0, 1, 4.5, 9, 99, 99, 99, 99)
  sd   <- 0.5
  nrep <- 4000
  
  set.seed(2024)
  q_R <- sobol_g2_qoi_mean_R(X, nrep = nrep, sd = sd, a = a)
  set.seed(2024)
  q_cpp <- sobol_g2_qoi_mean(X, nrep = nrep, sd = sd, a = a)
  
  expect_gt(cor(q_cpp, q_R), 0.999)
})

test_that("sobol_g2_with_covariate_noise_cpp reproduces R version under same seed", {
  set.seed(123)
  X <- cbind(
    matrix(runif(30 * 2), ncol = 2),
    runif(30, min = 1, max = 10)
  )
  a <- c(0, 1, 4.5, 9, 99, 99, 99, 99)
  
  set.seed(777)
  y_R <- sobol_g2_with_covariate_noise_R(X, a = a)
  set.seed(777)
  y_cpp <- sobol_g2_with_covariate_noise(X, a = a)
  
  expect_equal(y_cpp, y_R)
})

test_that("sobol_g2_qoi_covariate_mean_cpp reproduces R version under same seed", {
  set.seed(123)
  X <- cbind(
    matrix(runif(25 * 2), ncol = 2),
    runif(25, min = 1, max = 10)
  )
  a    <- c(0, 1, 4.5, 9, 99, 99, 99, 99)
  nrep <- 4000
  
  set.seed(999)
  q_R <- sobol_g2_qoi_covariate_mean_R(X, nrep = nrep, a = a)
  set.seed(999)
  q_cpp <- sobol_g2_qoi_covariate_mean(X, nrep = nrep, a = a)
  
  expect_gt(cor(q_cpp, q_R), 0.999)
})
