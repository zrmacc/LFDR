# Purpose: Data generation from univariate normal mixture.
# Updated: 180923

#' Simulated from Truncated Normal Distribution
#' 
#' Draws realizations from a truncated normal distribution.
#' 
#' @param n Observations.
#' @param L Lower truncation limit.
#' @param U Upper truncation limit.
#' @param m Location parameter.
#' @param s Scale parameter.
#' 
#' @return Vector of realizations from the truncated normal distribution. 
#' 
#' @importFrom stats runif pnorm qnorm
#' @export

rTN = function(n,L,U,m=0,s=1){
  # Alpha
  a = (L-m)/s;
  # Beta
  b = (U-m)/s;
  # Normalizing constant
  H = pnorm(b)-pnorm(a);
  # Uniforms
  u = runif(n=n);
  # Transform
  z = qnorm(pnorm(a)+H*u)*s+m;
  # Output
  return(z);
}

#' Simulated from Normal Mixture Model
#'
#' Draws realizations from a normal mixture model. First, the cluster
#' membership is drawn from a multinomial distribution, with mixture proportions
#' specified by \code{pi}. Conditional on cluster membership, the observation is
#' drawn from a multivariate normal distribution, with cluster-specific mean and
#' variance. The cluster means are provided using \code{m}, and the cluster
#' covariance matrices are provided using \code{v}.
#'
#' @param n Observations.
#' @param k Number of mixture components. Defaults to 1.
#' @param pi Mixture proportions. If omitted, components are assumed
#'   equi-probable.
#' @param m Vector of means, one for each component \eqn{k}.
#' @param v Vector of variances, one for each component \eqn{k}.
#'
#' @return Numeric matrix containing the mixture component and the observation.
#'
#' @importFrom foreach foreach '%do%'
#' @importFrom stats rbinom rmultinom rnorm 
#' @export
#' 
#' @examples
#' \dontrun{
#' # Two component mixture
#' Y = rNM(n=1e3,k=2,pi=c(0.95,0.05),m=c(0,2),v=c(1,1));
#' }

rNM = function(n,k=1,pi=NULL,m=0,v=1){
  ## Input checks
  # Mixture proportions
  if(is.null(pi)){pi = rep(1,k)/k} else {pi = pi/sum(pi)};
  if(length(pi)!=k){stop("Mixture proportions pi must have length k.")};
  # Means
  if(length(m)!=k){stop("One mean is expected for each mixture component k.")};
  # Variances
  if(length(v)!=k){stop("One variance is expected for each mixture component k.")};
  
  # Single component
  if(k==1){
    # Mixture component
    z = rep(1,n);
    # Draw response
    y = rnorm(n=n,mean=m,sd=sqrt(v));
  } else {
    # Draw mixture component
    z = rmultinom(n=n,size=1,prob=pi);
    aux = function(x){which(x==1)};
    z = apply(z,2,aux);
    # Loop over mixture components
    i = NULL;
    y = foreach(i=1:k,.combine=c) %do% {
      # Obs from current component 
      ni = sum(z==i);
      # Output only if component is non-empty
      if(ni>0){
        return(rnorm(n=ni,mean=m[i],sd=sqrt(v[i])));
      }
    };
  };
    
  # Output
  Out = cbind(z,y);
  colnames(Out) = c("Comp","Obs");
  rownames(Out) = seq(1:n);
  return(Out);
}