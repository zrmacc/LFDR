# Purpose: Fit an exponential family density.
# Updated: 2026-03-06.

#' Support
#'
#' Return an interval containing the input data.
#'
#' @param z Observations.
#' @param sig_digits Significant digits.
#'
#' @return Numeric vector containing a lower and upper bound.
Support <- function(z, sig_digits = 1) {
  # Data support (min and max).
  min_z <- min(z)
  max_z <- max(z)
  # Coarsen bounds to sig_digits.
  scale <- 10^(sig_digits)
  min_z <- floor(min_z * scale) / scale
  max_z <- ceiling(max_z * scale) / scale
  # Return named bounds.
  out <- c(min_z, max_z)
  names(out) <- c("L", "U")
  return(out)
}

########################
# Estimation function.
########################

#' Fit Exponential Family Density
#'
#' Fit an exponential family distribution of specified degree to a
#' sample of observations.
#'
#' @param z Observations.
#' @param n_bins Bins for discretizing.
#' @param degree Polynomial degree.
#' @param method Either "ns" for natural splines, or "poly" for orthogonal
#'   polynomials. Default is "poly".
#'
#' @return A list containing a function \code{f} for evaluating the fitted
#'   density, and the estimated regression coefficients \code{beta}.
#'
#' @export
FitEF <- function(
    z,
    n_bins = NULL,
    degree = 2,
    method = "poly"
) {
  ## Validate inputs.
  if (is.null(n_bins)) {
    n_bins <- ceiling(length(z) / 25)
  }
  # Validate method argument.
  choices <- c("ns", "poly")
  if (!(method %in% choices)) {
    stop("Select method from among: ns or poly.")
  }
  # Sample size.
  n_obs <- length(z)
  # Data range for binning.
  r <- Support(z)
  # Bin edges.
  bin <- seq(from = r[1], to = r[2], length.out = (n_bins + 1))
  # Width of each bin.
  bin_width <- (r[2] - r[1]) / n_bins
  # Bin midpoints for design matrix.
  mid_pts <- bin[1:(n_bins)] + (1 / 2) * diff(bin)
  # Quantile-based knot locations.
  loc <- stats::quantile(z, probs = seq(from = 0, to = degree) / degree)
  # Knots at data boundaries.
  loc_boundary <- c(loc[1], loc[degree + 1])
  # Interior knots (for splines).
  loc_interior <- c(loc[2:degree])
  # Assign observations to bins.
  zd <- cut(x = z, breaks = bin, include.lowest = TRUE)
  y <- as.numeric(table(zd))

  # Build regression basis (splines or polynomials).
  if (method == "ns") {
    # Design matrix with natural spline basis.
    x_mat <- cbind(1, splines::ns(mid_pts, knots = loc_interior, intercept = FALSE, Boundary.knots = loc_boundary))
    colnames(x_mat) <- paste0(rep("s", times = (degree + 1)), seq(from = 0, to = degree))
  } else if (method == "poly") {
    # Design matrix with orthogonal polynomials.
    x_mat <- stats::poly(mid_pts, degree = degree)
    # Keep poly coefs for later evaluation at new x.
    poly_coefs <- attributes(x_mat)$coefs
    # Add intercept column.
    x_mat <- cbind(1, x_mat)
    colnames(x_mat) <- paste0(rep("p", times = (degree + 1)), seq(from = 0, to = degree))
  }

  # Fit Poisson (count) regression.
  offset <- rep(n_obs * bin_width, times = n_bins)
  glm_fit <- stats::glm.fit(x = x_mat, y = y, offset = log(offset), family = stats::poisson(link = "log"))
  # Intercept and basis coefficients.
  beta <- stats::coef(glm_fit)
  b0 <- beta[1]
  b1 <- beta[2:(degree + 1)]
  # Closure that evaluates fitted density at new x.
  if (method == "ns") {
    f <- function(x) {
      p <- splines::ns(x, knots = loc_interior, intercept = FALSE, Boundary.knots = loc_boundary)
      out <- as.numeric(exp(b0 + c(p %*% b1)))
      return(out)
    }
  } else if (method == "poly") {
    f <- function(x) {
      out <- as.numeric(exp(b0 + c(stats::poly(x, degree = degree, coefs = poly_coefs) %*% b1)))
      return(out)
    }
  }
  # Return method, density function, and coefficients.
  out <- list("Method" = method, "f" = f, "beta" = beta)
  return(out)
}

########################
# Sampling function.
########################

#' Sample Fitted Exponential Family
#'
#' Sampling via the accept-reject algorithm. The proposal distribution
#' is taken as normal, with the option to specify the mean and variance.
#'
#' @param z Observations.
#' @param f Density to sample.
#' @param prop_mean Proposal mean.
#' @param prop_var Proposal variance.
#' @param n_samp Target number of samples to obtain.
#' @param parallel Ignored; retained for backward compatibility.
#' @param max_it Maximum number of accept-reject iterations.
#'
#' @return Vector of realizations from the mixture density.
#'
#' @export
SampleEF <- function(
    z,
    f,
    prop_mean = NULL,
    prop_var = NULL,
    n_samp = 1e3,
    parallel = FALSE,
    max_it = 25
) {
  # Data range for proposal support.
  support_r <- Support(z)
  # Default proposal mean and variance from data.
  if (is.null(prop_mean)) {
    prop_mean <- mean(z)
  }
  if (is.null(prop_var)) {
    prop_var <- stats::var(z)
  }
  # Normal proposal density function.
  prop_dens <- function(x) {
    stats::dnorm(x, mean = prop_mean, sd = sqrt(prop_var))
  }
  # Target density / proposal density (for accept-reject).
  lr <- function(x) {
    f(x) / prop_dens(x)
  }
  # Find maximum of likelihood ratio (upper bound for accept-reject).
  opt <- stats::optimize(f = lr, interval = support_r, maximum = TRUE)
  # Ceiling to 1 decimal for accept-reject bound M.
  m_lr <- ceiling(opt$objective * 10) / 10
  ## Accept-reject sampling loop.
  out_raw <- rep(NA_real_, n_samp)
  for (i in seq_len(n_samp)) {
    # Per-draw state.
    iter <- 0
    x <- NA_real_
    keep <- FALSE
    while (!keep) {
      # Draw from proposal distribution.
      z_prop <- stats::rnorm(n = 1, mean = prop_mean, sd = sqrt(prop_var))
      # Uniform variate for accept-reject comparison.
      u <- stats::runif(n = 1)
      # Accept if u < (1/M) * LR(z_prop).
      keep <- (u < (1 / m_lr) * lr(z_prop))
      # Store draw and exit loop when accepted.
      if (keep) {
        x <- z_prop
        break
      }
      # Increment iteration count when rejected.
      iter <- iter + 1
      if (iter >= max_it) {
        break
      }
    }
    out_raw[i] <- x
  }
  # Return accepted draws (drop NA from early exit).
  out <- as.numeric(out_raw[!is.na(out_raw)])
  return(out)
}
