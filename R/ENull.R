# Purpose: Estimation of empirical null distribution.
# Updated: 2026-03-06.

#' Estimate Empirical Null by Central Matching
#'
#' Estimates parameters of the empirical null distribution by central matching
#' within a "null neighborhood". Observations in the null neighborhood are
#' assumed to arise from the null component of the mixture distribution.
#'
#' @param z Observations.
#' @param f Fitted mixture density. See \code{\link{FitEF}}.
#' @param lower_bound Lower limit of null neighborhood.
#' @param upper_bound Upper limit of null neighborhood.
#'
#' @return Estimated mean and variance of the null component, and the proportion
#'   of observations arising from the null component.
#'

ENullCM <- function(z, f, lower_bound, upper_bound) {
  # Fine grid over null neighborhood.
  x <- seq(from = lower_bound, to = upper_bound, length.out = 1e3)
  # Intercept and quadratic terms.
  x_mat <- cbind(1, stats::poly(x, degree = 2, raw = TRUE, simple = TRUE))
  # Singular value decomposition of design matrix.
  svd_out <- svd(x_mat)
  u_mat <- svd_out$u
  v_mat <- svd_out$v
  d_vec <- svd_out$d
  # Log fitted mixture density on grid.
  y <- log(f(x))
  # Coefficients in orthogonalized (SVD) space.
  g0 <- t(u_mat) %*% y
  # Transform back to original parameterization.
  b0 <- as.numeric(v_mat %*% (g0 / d_vec))
  # Null variance from quadratic coefficient.
  v0 <- -1 / (2 * b0[3])
  # Null mean from linear coefficient.
  m0 <- v0 * b0[2]
  # Estimated proportion of null component.
  p0 <- exp(b0[1] + (1 / 2) * (m0^2 / v0 + log(2 * pi * v0)))
  # Return null mean, variance, and proportion.
  out <- list("m0" = m0, "v0" = v0, "p0" = p0)
  return(out)
}

#' Maximum Likelihood Estimation of the Empirical Null
#'
#' Estimates parameters of the empirical null distribution by fitting
#' a truncated normal distribution within a "null neighborhood". Observations
#' in the null neighborhood are assumed to arise from the null component
#' of the mixture distribution.
#'
#' @param z Observations.
#' @param lower_bound Lower limit of null neighborhood.
#' @param upper_bound Upper limit of null neighborhood.
#' @param max_it Maximum number of NR iterations.
#'
#' @return Estimated mean and variance of the null component, and the proportion
#'   of observations arising from the null component.
#'

ENullTN <- function(z, lower_bound, upper_bound, max_it = 10) {
  # Sample size.
  n_obs <- length(z)
  # Observations falling in null neighborhood.
  z0 <- z[z >= lower_bound & z <= upper_bound]
  n0 <- length(z0)
  # ML fit of truncated normal in null region.
  m_fit <- FitTN(
    z = z0, lower_bound = lower_bound, upper_bound = upper_bound,
    max_it = max_it, report = FALSE
  )
  # Fitted location (null mean).
  m0 <- m_fit[1]
  # Fitted scale (null SD).
  s0 <- m_fit[2]
  # Fraction of data in null region, renormalized by truncation constant.
  p0 <- n0 / (n_obs * H(
    loc_param = m0, scale_param = s0,
    lower_bound = lower_bound, upper_bound = upper_bound,
    order = 0, direction = "m"
  ))
  p0 <- min(p0, 1)
  # Return null mean, variance, and proportion.
  out <- list("m0" = m0, "v0" = s0^2, "p0" = p0)
  return(out)
}

#' Estimation of the Empirical Null
#'
#' Estimates the parameters of the empirical null distribution using
#' observations within a "null neighborhood". Observations in the null
#' neighborhood are assumed to arise from the null component of the mixture
#' distribution. The null neighborhood is specified using either a lower and
#' upper bound, or the central proportion of z scores supposed to arise from the
#' null.
#'
#' @param z Observations.
#' @param lower_bound Lower limit of null neighborhood.
#' @param upper_bound Upper limit of null neighborhood.
#' @param central_prop Central proportion of z scores belonging to null neighborhood.
#' @param method Either "CM" (central matching) or "ML" (maximum likelihood).
#' @param f Estimated mixture density if \code{method = "CM"}.
#'
#' @return List containing the estimated mean and variance of the null
#'   component, and the proportion of observations arising from the null
#'   component.
#'
#' @export
ENull <- function(
    z,
    lower_bound = NULL,
    upper_bound = NULL,
    central_prop = NULL,
    method = "ML",
    f = NULL) {
  ## Validate inputs.
  # Validate method and set default central_prop.
  if (!(method %in% c("CM", "ML"))) {
    stop("Select method from CM or ML.")
  }
  if (method == "CM" && is.null(f)) {
    stop("For central matching, provide an estimate of the mixture density.")
  }
  if (method == "CM" && is.null(central_prop)) {
    central_prop <- 0.5
  }
  if (method == "ML" && is.null(central_prop)) {
    central_prop <- 0.8
  }
  if (is.null(lower_bound)) {
    lower_bound <- as.numeric(stats::quantile(x = z, probs = (1 - central_prop) / 2))
  }
  if (is.null(upper_bound)) {
    upper_bound <- as.numeric(stats::quantile(x = z, probs = (1 + central_prop) / 2))
  }
  # Dispatch to CM or ML estimator.
  if (method == "CM") {
    theta <- ENullCM(
      z = z, f = f,
      lower_bound = lower_bound, upper_bound = upper_bound
    )
  } else {
    theta <- ENullTN(
      z = z,
      lower_bound = lower_bound, upper_bound = upper_bound
    )
  }
  # Return null mean, variance, and proportion.
  return(theta)
}
