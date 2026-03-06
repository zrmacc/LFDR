test_that("RTN generates within bounds", {
  set.seed(1)
  z <- RTN(n_obs = 1000, lower_bound = -2, upper_bound = 2, loc_param = 0, scale_param = 1)
  expect_true(all(z >= -2 & z <= 2))
  expect_length(z, 1000)
})

test_that("FitTN recovers parameters for truncated normal", {
  set.seed(2)
  z <- RTN(n_obs = 500, lower_bound = -1.5, upper_bound = 1.5, loc_param = 0.2, scale_param = 0.9)
  fit <- FitTN(z = z, lower_bound = -1.5, upper_bound = 1.5)
  expect_named(fit, c("m", "s"))
  expect_equal(unname(fit["m"]), 0.2, tolerance = 0.2)
  expect_equal(unname(fit["s"]), 0.9, tolerance = 0.2)
})

test_that("H renormalizing constant is positive and below 1", {
  h0 <- LFDR:::H(
    loc_param = 0, scale_param = 1,
    lower_bound = -2, upper_bound = 2,
    order = 0, direction = "m"
  )
  expect_true(h0 > 0 && h0 < 1)
})

test_that("QTN and UTN are consistent", {
  set.seed(3)
  z <- RTN(n_obs = 50, lower_bound = -1, upper_bound = 1, loc_param = 0, scale_param = 1)
  u <- LFDR:::UTN(
    z = z, loc_param = 0, log_scale = 0,
    lower_bound = -1, upper_bound = 1
  )
  expect_length(u, 2)
  expect_named(u, c("um", "ut"))
})

test_that("InfoTN returns 2x2 matrix", {
  set.seed(4)
  z <- RTN(n_obs = 30, lower_bound = -1, upper_bound = 1)
  info <- LFDR:::InfoTN(
    z = z, loc_param = 0, log_scale = 0,
    lower_bound = -1, upper_bound = 1
  )
  expect_equal(dim(info), c(2L, 2L))
  expect_true(all(is.finite(info)))
})
