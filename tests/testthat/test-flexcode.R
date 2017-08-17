context("Test FlexCoDE")

test_that("FlexCoDE fits previous results", {
  set.seed(42)

  n_obs <- 117
  n_x <- 3

  x_train <- matrix(rnorm(n_obs * n_x), n_obs, n_x)
  z_train <- matrix(rnorm(n_obs), n_obs, 1)

  x_validation <- matrix(rnorm(n_obs * n_x), n_obs, n_x)
  z_validation <- matrix(rnorm(n_obs), n_obs, 1)

  system <- "Fourier"
  n_basis <- 31

  obj <- fitFlexCoDE(x_train, z_train, x_validation, z_validation,
                     regressionFunction = regressionFunction.NN,
                     system = system, nIMax = n_basis)

  expect_equal_to_reference(obj, "cache/flexcode.rds")
})
