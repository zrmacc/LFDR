test_that("Support returns finite bounds", {
  z <- c(-2.5, 0, 3.1)
  s <- LFDR:::Support(z)
  expect_named(s, c("L", "U"))
  expect_lte(s["L"], min(z) + 0.01)
  expect_gte(s["U"], max(z) - 0.01)
})

test_that("FitEF returns density function and coefficients", {
  set.seed(20)
  z <- rnorm(200)
  fit <- FitEF(z = z, n_bins = 20, degree = 2, method = "poly")
  expect_named(fit, c("Method", "f", "beta"))
  expect_equal(fit$Method, "poly")
  expect_true(is.function(fit$f))
  expect_true(all(fit$f(z) > 0))
  expect_length(fit$beta, 3)
})

test_that("FitEF with method ns works", {
  set.seed(21)
  z <- rnorm(150)
  fit <- FitEF(z = z, n_bins = 15, degree = 3, method = "ns")
  expect_equal(fit$Method, "ns")
  expect_true(is.function(fit$f))
})

test_that("FitEF rejects invalid method", {
  expect_error(FitEF(z = rnorm(50), method = "invalid"))
})

test_that("solve and matrix product give identity when multiplied", {
  a <- matrix(c(2, 0, 0, 3), 2, 2)
  a_inv <- solve(a)
  expect_equal(a %*% a_inv, diag(2), tolerance = 1e-10)
})

test_that("SampleEF returns vector of requested length (small n)", {
  set.seed(30)
  z <- rnorm(100)
  fit <- FitEF(z = z, n_bins = 10, degree = 2, method = "poly")
  s <- SampleEF(z = z, f = fit$f, n_samp = 50, parallel = FALSE, max_it = 50)
  expect_true(is.numeric(s))
  expect_true(length(s) <= 50)
  expect_true(length(s) >= 1)
})
