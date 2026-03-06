test_that("FDR errors when f is missing", {
  z <- rnorm(50)
  expect_error(FDR(z = z, null_mean = 0, null_var = 1, null_prop = 1), "f.*required")
})

test_that("FDR returns lfdr and tFDR functions", {
  set.seed(50)
  z <- rnorm(200)
  fit <- FitEF(z = z, n_bins = 25, degree = 2, method = "poly")
  null_est <- ENull(z = z, lower_bound = -2, upper_bound = 2, method = "ML")
  out <- FDR(
    z = z,
    f = fit$f,
    n_bins = 25,
    null_mean = null_est$m0,
    null_var = null_est$v0,
    null_prop = null_est$p0
  )
  expect_named(out, c("lfdr", "tFDR"))
  expect_true(is.function(out$lfdr))
  expect_true(is.function(out$tFDR))
  lv <- out$lfdr(z)
  tv <- out$tFDR(z)
  expect_true(all(lv >= 0 & lv <= 1, na.rm = TRUE))
  expect_true(all(tv >= 0 & tv <= 1, na.rm = TRUE))
})

test_that("FDR works with default bins", {
  set.seed(51)
  z <- rnorm(150)
  fit <- FitEF(z = z, degree = 2, method = "poly")
  null_est <- ENull(z = z, lower_bound = -2, upper_bound = 2, method = "ML")
  out <- FDR(
    z = z, f = fit$f,
    null_mean = null_est$m0, null_var = null_est$v0, null_prop = null_est$p0
  )
  expect_named(out, c("lfdr", "tFDR"))
})

########################
# LFDRPipeline.
########################

test_that("LFDRPipeline returns fit, null, and fdr", {
  set.seed(60)
  z <- rnorm(200)
  out <- LFDRPipeline(
    z = z,
    n_bins = 30,
    degree = 2,
    fit_method = "poly",
    lower_bound = -2,
    upper_bound = 2,
    null_method = "ML"
  )
  expect_named(out, c("fit", "null", "fdr"))
  expect_named(out$fit, c("Method", "f", "beta"))
  expect_named(out$null, c("m0", "v0", "p0"))
  expect_named(out$fdr, c("lfdr", "tFDR"))
  expect_true(is.function(out$fit$f))
  expect_true(is.function(out$fdr$lfdr))
  expect_true(is.function(out$fdr$tFDR))
})

test_that("LFDRPipeline with default arguments runs", {
  set.seed(61)
  z <- rnorm(150)
  out <- LFDRPipeline(z = z)
  expect_named(out, c("fit", "null", "fdr"))
  expect_equal(length(out$null$m0), 1)
  expect_equal(length(out$null$v0), 1)
  expect_equal(length(out$null$p0), 1)
})

test_that("LFDRPipeline matches manual FitEF then ENull then FDR", {
  set.seed(62)
  z <- rnorm(180)
  n_bins <- 25
  degree <- 2
  lower_bound <- -2
  upper_bound <- 2
  pipe <- LFDRPipeline(
    z = z,
    n_bins = n_bins,
    degree = degree,
    fit_method = "poly",
    lower_bound = lower_bound,
    upper_bound = upper_bound,
    null_method = "ML"
  )
  fit_manual <- FitEF(z = z, n_bins = n_bins, degree = degree, method = "poly")
  null_manual <- ENull(
    z = z,
    lower_bound = lower_bound,
    upper_bound = upper_bound,
    method = "ML"
  )
  fdr_manual <- FDR(
    z = z,
    f = fit_manual$f,
    n_bins = n_bins,
    null_mean = null_manual$m0,
    null_var = null_manual$v0,
    null_prop = null_manual$p0
  )
  expect_equal(pipe$fit$Method, fit_manual$Method)
  expect_equal(pipe$null$m0, null_manual$m0)
  expect_equal(pipe$null$v0, null_manual$v0)
  expect_equal(pipe$null$p0, null_manual$p0)
  expect_equal(pipe$fdr$lfdr(z), fdr_manual$lfdr(z), tolerance = 1e-10)
  expect_equal(pipe$fdr$tFDR(z), fdr_manual$tFDR(z), tolerance = 1e-10)
})

test_that("LFDRPipeline works with central matching", {
  set.seed(63)
  z <- rnorm(200)
  out <- LFDRPipeline(
    z = z,
    n_bins = 40,
    degree = 3,
    fit_method = "poly",
    lower_bound = -2,
    upper_bound = 2,
    null_method = "CM"
  )
  expect_named(out$null, c("m0", "v0", "p0"))
  expect_true(is.finite(out$null$m0))
  expect_true(out$null$v0 > 0)
  expect_true(out$null$p0 > 0 && out$null$p0 <= 1)
})

########################
# FDRSummary.
########################

test_that("FDRSummary returns data frame with z, lfdr, tFDR", {
  set.seed(70)
  z <- rnorm(100)
  fit <- FitEF(z = z, n_bins = 20, degree = 2, method = "poly")
  null_est <- ENull(z = z, lower_bound = -2, upper_bound = 2, method = "ML")
  fdr <- FDR(
    z = z,
    f = fit$f,
    n_bins = 20,
    null_mean = null_est$m0,
    null_var = null_est$v0,
    null_prop = null_est$p0
  )
  tbl <- FDRSummary(z = z, fdr = fdr)
  expect_s3_class(tbl, "data.frame")
  expect_named(tbl, c("z", "lfdr", "tFDR"))
  expect_equal(nrow(tbl), length(z))
  expect_equal(tbl$z, z)
  expect_equal(tbl$lfdr, fdr$lfdr(z))
  expect_equal(tbl$tFDR, fdr$tFDR(z))
  expect_true(all(tbl$lfdr >= 0 & tbl$lfdr <= 1, na.rm = TRUE))
  expect_true(all(tbl$tFDR >= 0 & tbl$tFDR <= 1, na.rm = TRUE))
})

test_that("FDRSummary with cutoff adds discovery column", {
  set.seed(71)
  z <- rnorm(80)
  fit <- FitEF(z = z, n_bins = 15, degree = 2, method = "poly")
  null_est <- ENull(z = z, lower_bound = -2, upper_bound = 2, method = "ML")
  fdr <- FDR(
    z = z, f = fit$f,
    null_mean = null_est$m0, null_var = null_est$v0, null_prop = null_est$p0
  )
  tbl <- FDRSummary(z = z, fdr = fdr, cutoff = 0.2)
  expect_named(tbl, c("z", "lfdr", "tFDR", "discovery"))
  expect_equal(tbl$discovery, tbl$lfdr <= 0.2)
  expect_type(tbl$discovery, "logical")
})

test_that("FDRSummary discovery is TRUE where lfdr <= cutoff", {
  set.seed(72)
  z <- rnorm(50)
  fit <- FitEF(z = z, n_bins = 15, degree = 2, method = "poly")
  null_est <- ENull(z = z, lower_bound = -2, upper_bound = 2, method = "ML")
  fdr <- FDR(
    z = z, f = fit$f,
    null_mean = null_est$m0, null_var = null_est$v0, null_prop = null_est$p0
  )
  tbl <- FDRSummary(z = z, fdr = fdr, cutoff = 0.5)
  expect_equal(tbl$discovery, tbl$lfdr <= 0.5)
})

test_that("FDRSummary errors when fdr is invalid", {
  z <- rnorm(10)
  expect_error(
    FDRSummary(z = z, fdr = list(lfdr = identity)),
    "fdr must be the output of FDR"
  )
  expect_error(
    FDRSummary(z = z, fdr = list(tFDR = identity)),
    "fdr must be the output of FDR"
  )
  expect_error(
    FDRSummary(z = z, fdr = "not a list"),
    "fdr must be the output of FDR"
  )
})
