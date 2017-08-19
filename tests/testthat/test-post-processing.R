context("Test post-processing functions")


test_that("normalize_density deflates densities", {
  set.seed(42)

  density <- runif(100)
  bin_size <- 1

  normalize <- normalize_density(bin_size, density)

  expect_true(all(normalize >= 0.0))
  expect_equal(bin_size * sum(normalize), 1, tolerance = 1e-3)
  expect_equal(which.max(normalize), which.max(density))
})

test_that("normalize_density inflates densities", {
  set.seed(43)

  density <- runif(100)
  bin_size <- 1 / 100

  normalize <- normalize_density(bin_size, density)

  expect_true(all(normalize >= 0.0))
  expect_equal(bin_size * sum(normalize), 1, tolerance = 1e-3)
  expect_equal(which.max(normalize), which.max(density))
})
