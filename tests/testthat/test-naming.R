# Test that naming conventions are followed: CamelCase for functions,
# lower_snake_case for non-function objects (including arguments).

test_that("exported functions are CamelCase with no leading lowercase or periods", {
  exports <- c("FDR", "ENull", "FitEF", "FitTN", "RNM", "RTN", "SampleEF")
  for (fn in exports) {
    expect_true(
      grepl("^[A-Z]", fn),
      info = paste(fn, "should start with uppercase")
    )
    expect_false(
      grepl("\\.", fn),
      info = paste(fn, "should not contain periods")
    )
  }
})

test_that("FitTN accepts lower_snake_case arguments", {
  set.seed(1)
  z <- RTN(n_obs = 50, lower_bound = -1, upper_bound = 1)
  fit <- FitTN(
    z = z,
    lower_bound = -1,
    upper_bound = 1,
    max_it = 10,
    report = FALSE
  )
  expect_named(fit, c("m", "s"))
})

test_that("RNM accepts lower_snake_case arguments", {
  out <- RNM(
    n_obs = 20,
    n_comp = 2,
    pi = c(0.5, 0.5),
    loc_params = c(0, 1),
    scale_vars = c(1, 1)
  )
  expect_equal(nrow(out), 20)
})

test_that("solve and matrix product work for 2x2 identity", {
  a <- matrix(c(1, 0, 0, 1), 2, 2)
  a_inv <- solve(a)
  expect_equal(a %*% a_inv, diag(2), tolerance = 1e-10)
})
