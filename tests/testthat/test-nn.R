context("Test NN regression")

test_that("NN fits previous results", {
  set.seed(42)

  n_obs <- 110
  n_x <- 7
  n_z <- 13

  x <- matrix(rnorm(n_obs * n_x), n_obs, n_x)
  z <- matrix(rnorm(n_obs * n_z), n_obs, n_z)

  obj <- regressionFunction.NN(x, z)

  expect_equal_to_reference(obj, "cache/nn-regression.rds")
})

test_that("NN predict fits previous results", {
  set.seed(42)

  n_train <- 110
  n_pred <- 103
  n_x <- 7
  n_z <- 13

  x <- matrix(rnorm(n_train * n_x), n_train, n_x)
  x_pred <- matrix(rnorm(n_pred * n_x), n_pred, n_x)
  z <- matrix(rnorm(n_train * n_z), n_train, n_z)

  reg <- regressionFunction.NN(x, z)
  obj <- predict(reg, x_pred)

  expect_equal_to_reference(obj, "cache/nn-predict.rds")
})
