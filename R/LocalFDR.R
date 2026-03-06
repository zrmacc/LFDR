# Purpose: Calculate false discovery rates.
# Updated: 2026-03-06.

#' Local False Discovery Rate
#'
#' Estimates the local false discovery rate from a collection of z scores.
#'
#' @param z Observations.
#' @param f Mixture density function (required).
#' @param n_bins Bins.
#' @param null_mean Mean of empirical null.
#' @param null_var Variance of empirical null.
#' @param null_prop Null proportion.
#'
#' @return List including two functions, lfdr for calculating the local false
#'   discovery rate, and tFDR for calculating the tail false discovery rate.
#'
#' @export

FDR <- function(
    z,
    f = NULL,
    n_bins = NULL,
    null_mean = 0,
    null_var = 1,
    null_prop = 1) {
  if (is.null(f)) {
    stop("Argument f (mixture density function) is required.")
  }
  # Default number of bins if not specified.
  if (is.null(n_bins)) {
    n_bins <- ceiling(length(z) / 25)
  }
  # Sample size.
  n_obs <- length(z)
  # Data support from Support().
  r <- Support(z)
  range_width <- as.numeric(r[2] - r[1])
  # Bin edges.
  bin <- seq(from = r[1], to = r[2], length.out = (n_bins + 1))
  # Width of each bin.
  bin_width <- (range_width / n_bins)
  # Bin midpoints.
  mid_pts <- bin[1:(n_bins)] + (1 / 2) * diff(bin)
  # Matrix of bin intervals (left, right) for reference.
  a_mat <- cbind(bin[1:n_bins], bin[2:(n_bins + 1)])
  # Assign observations to bins.
  zd <- cut(x = z, breaks = bin, include.lowest = TRUE, right = FALSE)
  # Null SD (reused for CDF and tail integrals).
  null_sd <- sqrt(null_var)
  # Null CDF mass in each bin.
  q <- stats::pnorm(bin, mean = null_mean, sd = null_sd)
  q <- q[2:(n_bins + 1)] - q[1:n_bins]
  # Expected counts under null in each bin.
  expected_count <- n_obs * (null_prop * q)
  # Fitted mixture density times bin width (expected under mixture).
  y <- n_obs * bin_width * f(mid_pts)
  # Raw local FDR (expected/observed) per bin, capped at 1.
  aux <- function(x) {
    min(x, 1)
  }
  lfdr_vec <- vapply(expected_count / y, FUN = aux, FUN.VALUE = numeric(1))
  # Data frame of bin midpoints and local FDR values.
  curve_df <- data.frame("zd" = levels(zd), "mp" = mid_pts, "lfdr" = lfdr_vec)
  # Interpolating function for local FDR.
  lfdr <- stats::approxfun(
    x = curve_df$mp, y = curve_df$lfdr,
    method = "linear", f = 0, rule = 2
  )
  # Slightly extended range for tail integrals.
  lower_int <- min(bin) - 0.05 * range_width
  upper_int <- max(bin) + 0.05 * range_width
  # Tail FDR at each bin midpoint.
  tfdr <- numeric(n_bins)
  for (i in seq_len(n_bins)) {
    # Current bin midpoint.
    x <- mid_pts[i]
    if (x < 0) {
      s0 <- stats::pnorm(q = x, mean = null_mean, sd = null_sd)
      s_val <- stats::integrate(f = f, lower = lower_int, upper = x)$value
    } else {
      s0 <- stats::pnorm(q = x, mean = null_mean, sd = null_sd, lower.tail = FALSE)
      s_val <- stats::integrate(f = f, lower = x, upper = upper_int)$value
    }
    # Tail FDR at this point (capped at 1).
    tfdr[i] <- min(null_prop * s0 / s_val, 1)
  }
  # Interpolating function for tail FDR.
  t_fdr <- stats::approxfun(
    x = curve_df$mp, y = tfdr,
    method = "linear", f = 0, rule = 2
  )
  # Return local FDR and tail FDR functions.
  out <- list("lfdr" = lfdr, "tFDR" = t_fdr)
  return(out)
}

########################
# End-to-end pipeline.
########################

#' End-to-End Local FDR Pipeline
#'
#' Runs the full workflow: fit mixture density (\code{\link{FitEF}}), estimate
#' empirical null (\code{\link{ENull}}), then compute local and tail FDR
#' (\code{\link{FDR}}). Returns all intermediate results for inspection or
#' further use.
#'
#' @param z Observations (z-scores).
#' @param n_bins Number of bins for density and FDR estimation. Defaults to
#'   \code{ceiling(length(z) / 25)} if \code{NULL}.
#' @param degree Polynomial degree for mixture fit.
#' @param fit_method Basis for mixture fit: \code{"ns"} (natural splines) or
#'   \code{"poly"} (orthogonal polynomials).
#' @param lower_bound Lower limit of null neighborhood. If \code{NULL}, set from
#'   \code{central_prop} quantiles.
#' @param upper_bound Upper limit of null neighborhood. If \code{NULL}, set from
#'   \code{central_prop} quantiles.
#' @param central_prop Central proportion of z-scores defining the null
#'   neighborhood when bounds are \code{NULL}. Defaults to 0.5 for central
#'   matching and 0.8 for ML.
#' @param null_method Method for empirical null: \code{"CM"} (central matching)
#'   or \code{"ML"} (maximum likelihood). Central matching requires the fitted
#'   density and is used automatically in the pipeline.
#'
#' @return A list with components \code{fit} (output of \code{\link{FitEF}}),
#'   \code{null} (output of \code{\link{ENull}}), and \code{fdr} (output of
#'   \code{\link{FDR}}).
#'
#' @export
LFDRPipeline <- function(
    z,
    n_bins = NULL,
    degree = 2,
    fit_method = "poly",
    lower_bound = NULL,
    upper_bound = NULL,
    central_prop = NULL,
    null_method = "ML") {
  # Fit mixture density.
  fit <- FitEF(z = z, n_bins = n_bins, degree = degree, method = fit_method)
  # Estimate empirical null (uses fit$f when null_method is "CM").
  null <- ENull(
    z = z,
    lower_bound = lower_bound,
    upper_bound = upper_bound,
    central_prop = central_prop,
    method = null_method,
    f = fit$f
  )
  # Compute local and tail FDR.
  fdr <- FDR(
    z = z,
    f = fit$f,
    n_bins = n_bins,
    null_mean = null$m0,
    null_var = null$v0,
    null_prop = null$p0
  )
  out <- list("fit" = fit, "null" = null, "fdr" = fdr)
  return(out)
}

########################
# Summary table.
########################

#' FDR Summary Table
#'
#' Builds a data frame of observations with their local and tail FDR values, and
#' an optional discovery flag at a given lfdr cutoff.
#'
#' @param z Observations (z-scores).
#' @param fdr Output of \code{\link{FDR}}, a list with components \code{lfdr}
#'   and \code{tFDR} (functions).
#' @param cutoff Optional lfdr cutoff. If provided, a logical column
#'   \code{discovery} is added: \code{TRUE} where \code{lfdr(z) <= cutoff}.
#'
#' @return A data frame with columns \code{z}, \code{lfdr}, \code{tFDR}, and
#'   \code{discovery} (if \code{cutoff} is not \code{NULL}).
#'
#' @export
FDRSummary <- function(z, fdr, cutoff = NULL) {
  if (!is.list(fdr) || !all(c("lfdr", "tFDR") %in% names(fdr))) {
    stop("fdr must be the output of FDR(), with components lfdr and tFDR.")
  }
  lfdr_vals <- fdr$lfdr(z)
  tfdr_vals <- fdr$tFDR(z)
  out <- data.frame(
    "z" = z,
    "lfdr" = lfdr_vals,
    "tFDR" = tfdr_vals,
    stringsAsFactors = FALSE
  )
  if (!is.null(cutoff)) {
    out$discovery <- lfdr_vals <= cutoff
  }
  return(out)
}
