test_that("ENull ML returns m0, v0, p0", {
  set.seed(40)
  z <- rnorm(300)
  theta <- ENull(z = z, lower_bound = -2, upper_bound = 2, method = "ML")
  expect_named(theta, c("m0", "v0", "p0"))
  expect_equal(length(theta$m0), 1L)
  expect_equal(length(theta$v0), 1L)
  expect_equal(length(theta$p0), 1L)
  expect_gt(theta$v0, 0)
  expect_gte(theta$p0, 0)
  expect_lte(theta$p0, 1)
})

test_that("ENull CM requires f", {
  expect_error(ENull(z = rnorm(50), method = "CM"))
})

test_that("ENull CM works when f and lower_bound, upper_bound provided", {
  set.seed(41)
  z <- rnorm(200)
  fit <- FitEF(z = z, n_bins = 20, degree = 2, method = "poly")
  theta <- ENull(
    z = z, lower_bound = -1.5, upper_bound = 1.5,
    method = "CM", f = fit$f
  )
  expect_named(theta, c("m0", "v0", "p0"))
  expect_true(is.finite(theta$m0))
  expect_true(theta$v0 > 0)
})

test_that("ENull rejects invalid method", {
  expect_error(ENull(z = rnorm(50), method = "invalid"))
})
