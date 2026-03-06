# Purpose: Data generation from univariate normal mixture.
# Updated: 2026-03-06.

#' Simulate from Truncated Normal Distribution
#'
#' Draws realizations from a truncated normal distribution.
#'
#' @param n_obs Number of observations.
#' @param lower_bound Lower truncation limit.
#' @param upper_bound Upper truncation limit.
#' @param loc_param Location parameter.
#' @param scale_param Scale parameter.
#'
#' @return Vector of realizations from the truncated normal distribution.
#'
#' @export

RTN <- function(n_obs, lower_bound, upper_bound, loc_param = 0, scale_param = 1) {
  # Standardized lower limit.
  a <- (lower_bound - loc_param) / scale_param
  # Standardized upper limit.
  b <- (upper_bound - loc_param) / scale_param
  # CDF difference (truncation constant).
  h_val <- stats::pnorm(b) - stats::pnorm(a)
  # Uniform draws for inverse CDF transform.
  u <- stats::runif(n = n_obs)
  # Inverse CDF transform to truncated normal.
  z <- stats::qnorm(stats::pnorm(a) + h_val * u) * scale_param + loc_param
  # Return draws.
  return(z)
}

#' Simulate from Normal Mixture Model
#'
#' Draws realizations from a normal mixture model. First, the cluster
#' membership is drawn from a multinomial distribution, with mixture proportions
#' specified by \code{pi}. Conditional on cluster membership, the observation is
#' drawn from a normal distribution, with cluster-specific mean and
#' variance. The cluster means are provided using \code{loc_params}, and the
#' cluster variances are provided using \code{scale_vars}.
#'
#' @param n_obs Number of observations.
#' @param n_comp Number of mixture components. Defaults to 1.
#' @param pi Mixture proportions. If omitted, components are assumed
#'   equi-probable.
#' @param loc_params Vector of means, one for each component.
#' @param scale_vars Vector of variances, one for each component.
#'
#' @return Numeric matrix containing the mixture component and the observation.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Two component mixture
#' Y <- RNM(n_obs = 1e3, n_comp = 2, pi = c(0.95, 0.05),
#'          loc_params = c(0, 2), scale_vars = c(1, 1))
#' }
RNM <- function(
    n_obs,
    n_comp = 1,
    pi = NULL,
    loc_params = 0,
    scale_vars = 1) {
  ## Validate inputs.
  # Normalize and validate pi.
  if (is.null(pi)) {
    pi <- rep(1, n_comp) / n_comp
  } else {
    pi <- pi / sum(pi)
  }
  if (length(pi) != n_comp) {
    stop("Mixture proportions pi must have length n_comp.")
  }
  # Validate loc_params length.
  if (length(loc_params) != n_comp) {
    stop("One mean is expected for each mixture component.")
  }
  # Validate scale_vars length.
  if (length(scale_vars) != n_comp) {
    stop("One variance is expected for each mixture component.")
  }

  # k = 1: no mixture, just one normal.
  if (n_comp == 1) {
    # All observations from component 1.
    z <- rep(1, n_obs)
    # Draw from single normal.
    y <- stats::rnorm(n = n_obs, mean = loc_params, sd = sqrt(scale_vars))
  } else {
    # Multinomial component labels.
    z_mat <- stats::rmultinom(n = n_obs, size = 1, prob = pi)
    aux <- function(x) {
      which(x == 1)
    }
    z <- apply(z_mat, 2, aux)
    # Draw y[i] from component z[i] so observation order is preserved.
    y <- numeric(n_obs)
    for (i in seq_len(n_obs)) {
      y[i] <- stats::rnorm(
        n = 1,
        mean = loc_params[z[i]],
        sd = sqrt(scale_vars[z[i]])
      )
    }
  }

  # Return matrix with component and observation columns.
  out <- cbind(z, y)
  colnames(out) <- c("Comp", "Obs")
  rownames(out) <- seq_len(n_obs)
  return(out)
}
