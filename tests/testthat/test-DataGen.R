test_that("RNM single component returns correct structure", {
  set.seed(10)
  out <- RNM(n_obs = 100, n_comp = 1, loc_params = 0, scale_vars = 1)
  expect_equal(dim(out), c(100, 2))
  expect_equal(colnames(out), c("Comp", "Obs"))
  expect_true(all(out[, "Comp"] == 1))
})

test_that("RNM two components preserves observation order", {
  set.seed(11)
  out <- RNM(
    n_obs = 200, n_comp = 2, pi = c(0.5, 0.5),
    loc_params = c(-1, 1), scale_vars = c(1, 1)
  )
  expect_equal(dim(out), c(200, 2))
  comp <- out[, "Comp"]
  obs <- out[, "Obs"]
  for (i in seq_len(200)) {
    expect_true(comp[i] %in% c(1, 2))
  }
  expect_true(mean(obs[comp == 1]) < 0)
  expect_true(mean(obs[comp == 2]) > 0)
})

test_that("RNM rejects invalid pi length", {
  expect_error(RNM(
    n_obs = 10, n_comp = 2, pi = c(1, 1, 1),
    loc_params = c(0, 0), scale_vars = c(1, 1)
  ))
})

test_that("RNM rejects invalid loc_params length", {
  expect_error(RNM(
    n_obs = 10, n_comp = 2, pi = c(0.5, 0.5),
    loc_params = c(0), scale_vars = c(1, 1)
  ))
})

test_that("RNM rejects invalid scale_vars length", {
  expect_error(RNM(
    n_obs = 10, n_comp = 2, pi = c(0.5, 0.5),
    loc_params = c(0, 0), scale_vars = c(1)
  ))
})
