# Purpose: Estimated parameters for truncated normal distribution.
# Updated: 2026-03-06.

########################
# Renormalizing constant.
########################

# Formulas verified numerically.

#' Truncated Normal Renormalizing Constant
#'
#' Calculates the renormalizing constant of the truncated
#' normal distribution, and its derivatives.
#'
#' @param loc_param Location parameter.
#' @param scale_param Scale parameter.
#' @param lower_bound Lower truncation limit.
#' @param upper_bound Upper truncation limit.
#' @param order Order.
#' @param direction Direction.
#'

H <- function(loc_param, scale_param, lower_bound, upper_bound, order = 0, direction) {
  # Standardized lower and upper limits (alpha and beta).
  a <- (lower_bound - loc_param) / scale_param
  b <- (upper_bound - loc_param) / scale_param
  # Initialize output.
  out <- 0
  # Branch on derivative order.
  if (order == 0) {
    # Zeroth derivative (CDF difference).
    out <- stats::pnorm(b) - stats::pnorm(a)
  } else if (order == 1) {
    ## First derivatives.
    # Derivative with respect to location.
    if (direction == "m") {
      out <- -(1 / scale_param) * (stats::dnorm(b) - stats::dnorm(a))
    }
    # Derivative with respect to scale.
    if (direction == "s") {
      out <- -(1 / scale_param) * (b * stats::dnorm(b) - a * stats::dnorm(a))
    }
  } else if (order == 2) {
    ## Second derivatives.
    # Second derivative with respect to location.
    if (direction == "m") {
      out <- -(1 / scale_param)^2 * (b * stats::dnorm(b) - a * stats::dnorm(a))
    }
    # Second derivative with respect to scale.
    if (direction == "s") {
      out <- -(1 / scale_param)^2 * (b * (b^2 - 2) * stats::dnorm(b) - a * (a^2 - 2) * stats::dnorm(a))
    }
    # Mixed partial (location and scale).
    if (direction == "ms") {
      out <- -(1 / scale_param)^2 * ((b^2 - 1) * stats::dnorm(b) - (a^2 - 1) * stats::dnorm(a))
    }
  }
  # Return result.
  return(out)
}
########################
# Log likelihood.
########################

#' Truncated Normal Log Likelihood
#'
#' Evaluates the sample log likelihood of the truncated
#' normal distribution.
#'
#' @param z Observations.
#' @param loc_param Location parameter.
#' @param log_scale Log of scale parameter.
#' @param lower_bound Lower truncation limit.
#' @param upper_bound Upper truncation limit.
#'

QTN <- function(z, loc_param, log_scale, lower_bound, upper_bound) {
  # Recover scale parameter from log-scale.
  scale_param <- exp(log_scale)
  # Number of observations.
  n_obs <- length(z)
  # Normalizing constant (CDF over truncation interval).
  h0 <- H(
    loc_param = loc_param, scale_param = scale_param,
    lower_bound = lower_bound, upper_bound = upper_bound,
    order = 0, direction = "m"
  )
  # Log-likelihood value.
  out <- -n_obs * log(scale_param) - n_obs * log(h0) -
    sum((z - loc_param)^2) / (2 * scale_param^2)
  # Return result.
  return(out)
}
########################
# Score equations.
########################

# Formulas verified numerically.

#' Truncated Normal Score Equations
#'
#' Evaluates the gradient of the sample log likelihood.
#'
#' @param z Observations.
#' @param loc_param Location parameter.
#' @param log_scale Log of scale parameter.
#' @param lower_bound Lower truncation limit.
#' @param upper_bound Upper truncation limit.
#'

UTN <- function(z, loc_param, log_scale, lower_bound, upper_bound) {
  # Number of observations.
  n_obs <- length(z)
  # Recover scale parameter from log-scale.
  scale_param <- exp(log_scale)
  # Normalizing constant (CDF over truncation interval).
  h0 <- H(
    loc_param = loc_param, scale_param = scale_param,
    lower_bound = lower_bound, upper_bound = upper_bound,
    order = 0, direction = "m"
  )
  # First derivative of H with respect to location.
  h_m <- H(
    loc_param = loc_param, scale_param = scale_param,
    lower_bound = lower_bound, upper_bound = upper_bound,
    order = 1, direction = "m"
  )
  # First derivative of H with respect to scale.
  h_s <- H(
    loc_param = loc_param, scale_param = scale_param,
    lower_bound = lower_bound, upper_bound = upper_bound,
    order = 1, direction = "s"
  )
  # Score (gradient) component for location.
  u_m <- -(n_obs / h0) * h_m + sum((z - loc_param)) / (scale_param^2)
  # Score (gradient) component for log-scale.
  u_t <- (-(n_obs / scale_param) - (n_obs / h0) * h_s +
    sum((z - loc_param)^2) / (scale_param^3)) * scale_param
  # Name the score vector components.
  out <- c(u_m, u_t)
  names(out) <- c("um", "ut")
  # Return result.
  return(out)
}

########################
# Information matrix.
########################

# Formulas verified numerically.

#' Truncated Normal Observed Information
#'
#' Evaluates the observed information of the sample log likelihood.
#'
#' @param z Observations.
#' @param loc_param Location parameter.
#' @param log_scale Log of scale parameter.
#' @param lower_bound Lower truncation limit.
#' @param upper_bound Upper truncation limit.
#'

InfoTN <- function(z, loc_param, log_scale, lower_bound, upper_bound) {
  # Number of observations.
  n_obs <- length(z)
  # Recover scale parameter from log-scale.
  scale_param <- exp(log_scale)
  # Normalizing constant (CDF over truncation interval).
  h0 <- H(
    loc_param = loc_param, scale_param = scale_param,
    lower_bound = lower_bound, upper_bound = upper_bound,
    order = 0, direction = "m"
  )
  # First and second derivatives of H with respect to location.
  h_m <- H(
    loc_param = loc_param, scale_param = scale_param,
    lower_bound = lower_bound, upper_bound = upper_bound,
    order = 1, direction = "m"
  )
  h_mm <- H(
    loc_param = loc_param, scale_param = scale_param,
    lower_bound = lower_bound, upper_bound = upper_bound,
    order = 2, direction = "m"
  )
  # First and second derivatives of H with respect to scale.
  h_s <- H(
    loc_param = loc_param, scale_param = scale_param,
    lower_bound = lower_bound, upper_bound = upper_bound,
    order = 1, direction = "s"
  )
  h_ss <- H(
    loc_param = loc_param, scale_param = scale_param,
    lower_bound = lower_bound, upper_bound = upper_bound,
    order = 2, direction = "s"
  )
  # Mixed partial (location and scale).
  h_ms <- H(
    loc_param = loc_param, scale_param = scale_param,
    lower_bound = lower_bound, upper_bound = upper_bound,
    order = 2, direction = "ms"
  )
  # Observed information for location parameter.
  i_mm <- -n_obs / (h0^2) * (h_m^2) + (n_obs / h0) * h_mm + n_obs / (scale_param^2)
  # Cross-term (location and log-scale).
  i_ms <- -n_obs / (h0^2) * (h_m * h_s) + n_obs / h0 * h_ms +
    2 * sum((z - loc_param)) / (scale_param^3)
  i_mt <- i_ms * scale_param
  # Reuse score component for log-scale.
  u_t <- UTN(
    z = z, loc_param = loc_param, log_scale = log_scale,
    lower_bound = lower_bound, upper_bound = upper_bound
  )[2]
  # Observed information for log-scale parameter.
  h_ss_val <- n_obs / (scale_param^2) + n_obs / (h0^2) * (h_s^2) -
    n_obs / h0 * (h_ss) - 3 * sum((z - loc_param)^2) / (scale_param^4)
  i_tt <- -1 * (h_ss_val * scale_param^2 + u_t)
  # Return 2 x 2 information matrix.
  out <- matrix(c(i_mm, i_mt, i_mt, i_tt), nrow = 2)
  rownames(out) <- colnames(out) <- c("m", "t")
  return(out)
}

########################
# Fit truncated normal.
########################

#' Estimate Mean and Scale of Truncated Normal
#'
#' Estimates the location and scale of a truncated normal distribution.
#'
#' @param z Observations.
#' @param lower_bound Lower truncation limit.
#' @param upper_bound Upper truncation limit.
#' @param max_it Maximum iterations.
#' @param eps Tolerance.
#' @param report Report fitting progress?
#'
#' @return Vector containing the estimated location and scale of the truncated
#'   normal distribution.
#'
#' @export

FitTN <- function(
    z,
    lower_bound,
    upper_bound,
    max_it = 25,
    eps = 1e-8,
    report = FALSE) {
  # One Newton-Raphson step.
  update_fn <- function(theta) {
    # Current location and log-scale.
    m0 <- theta$m
    t0 <- theta$t
    # Log-likelihood at current parameters.
    q0 <- QTN(
      z = z, loc_param = m0, log_scale = t0,
      lower_bound = lower_bound, upper_bound = upper_bound
    )
    # Score vector at current parameters.
    u0 <- UTN(
      z = z, loc_param = m0, log_scale = t0,
      lower_bound = lower_bound, upper_bound = upper_bound
    )
    # Observed information matrix at current parameters.
    j0 <- InfoTN(
      z = z, loc_param = m0, log_scale = t0,
      lower_bound = lower_bound, upper_bound = upper_bound
    )
    # Newton-Raphson update (parameter proposal).
    prop <- c(m0, t0) + as.numeric(solve(j0) %*% u0)
    m1 <- prop[1]
    t1 <- prop[2]
    # Log-likelihood at proposed parameters.
    q1 <- QTN(
      z = z, loc_param = m1, log_scale = t1,
      lower_bound = lower_bound, upper_bound = upper_bound
    )
    # Log-likelihood increment (for convergence check).
    d <- q1 - q0
    # Return updated parameters and increment.
    out <- list("m" = m1, "t" = t1, "d" = d)
    return(out)
  }

  # Start at location 0, log-scale 0.
  theta0 <- list("m" = 0, "t" = 0)
  ## Maximization loop.
  for (i in seq_len(max_it)) {
    # Perform one Newton-Raphson step.
    theta1 <- update_fn(theta0)
    # Accept proposal if it increases log-likelihood.
    if (theta1$d > 0) {
      theta0 <- theta1
      if (report) {
        cat("Objective increment: ", signif(theta1$d, digits = 3), "\n")
      }
    }
    # Stop when improvement is negligible.
    if (theta1$d < eps) {
      break
    }
  }

  # Return location and scale (exp(log-scale)).
  out <- c(theta0$m, exp(theta0$t))
  names(out) <- c("m", "s")
  return(out)
}
