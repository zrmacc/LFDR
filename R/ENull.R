# Purpose: Estimation of empirical null distribution
# Updated: 180923

#' Estimate Empirical Null by Central Matching
#'
#' Estimates parameters of the empirical null distribution by central matching
#' within a "null neighborhood". Observations in the null neighborhood are
#' assumed to arise from the null component of the mixture distribution.
#' 
#' @param z Observations. 
#' @param f Fitted mixture density. See \code{\link{Fit.EF}}.
#' @param L Lower limit of null neighborhood.
#' @param U Upper limit of null neighborhood.
#' 
#' @return Estimated mean and variance of the null component, and the proportion
#'   of observations arising from the null component.
#' 
#' @importFrom stats coef lm.fit poly

eNull.CM = function(z,f,L,U){
  # Evaluation grid
  x = seq(from=L,to=U,length.out=1e3);
  # Design matrix
  X = cbind(1,poly(x,degree=2,raw=T,simple=T));
  # SVD
  S = svd(X);
  U = S$u;
  V = S$v;
  d = S$d;
  # Log mixture density
  y = log(f(x));
  # Regression coefficient w.r.t. orthogonal design
  g0 = matIP(U,y);
  # Convert to original scale
  b0 = as.numeric(MMP(V,g0/d));
  # Variance
  v0 = -1/(2*b0[3]);
  # Mean
  m0 = v0*b0[2];
  # Proportion null
  p0 = exp(b0[1]+(1/2)*(m0^2/v0+log(2*pi*v0)));
  # Output
  Out = list("m0"=m0,"v0"=v0,"p0"=p0);
  return(Out);
}

#' Maximum Likelihood Estimation of the Empirical Null
#'
#' Estimates parameters of the empirical null distribution by fitting
#' a truncated normal distribution within a "null neighborhood". Observations
#' in the null neighborhood are assumed to arise from the null component
#' of the mixture distribution. 
#' 
#' @param z Observations.
#' @param L Lower limit of null neighborhood.
#' @param U Upper limit of null neighborhood.
#' @param maxit Maximum number of NR iterations. 
#' 
#' @return Estimated mean and variance of the null component, and the proportion
#'   of observations arising from the null component.
#' 
#' @importFrom stats pnorm

eNull.TN = function(z,L,U,maxit=10){
  # Total observations
  n = length(z);
  # Obs in null neighborhood
  z0 = z[z>=L & z<=U];
  n0 = length(z0);
  # Fit truncated normal to obs in null neighborhood
  M = fitTN(z=z0,L=L,U=U,maxit=maxit,report=F); 
  # Location
  m0 = M[1];
  # Scale
  s0 = M[2];
  # Null proportion
  p0 = n0/(n*H(m=m0,s=s0,L=L,U=U,k=0));
  p0 = min(p0,1);
  # Output
  Out = list("m0"=m0,"v0"=s0^2,"p0"=p0);
  return(Out);
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
#' @param L Lower limit of null neighborhood.
#' @param U Upper limit of null neighborhood.
#' @param p Central proportion of z scores belonging to null neighborhood.
#' @param method Either 1. CM central matching, or 2. ML maximum likelihood
#' @param f Estimated mixture density if \code{method="CM"}.
#' 
#' @return List containing the estimated mean and variance of the null
#'   component, and the proportion of observations arising from the null
#'   component.
#' 
#' @importFrom stats quantile
#' @export 

eNull = function(z,L=NULL,U=NULL,p=NULL,method="ML",f=NULL){
  ## Input check
  # Method
  if(!(method %in% c("CM","ML"))){stop("Select method from CM or ML.")};
  if(method=="CM"&is.null(f)){stop("For central matching, provide an estimate of the mixture density.")};
  if(method=="CM"&is.null(p)){p=0.5};
  if(method=="ML"&is.null(p)){p=0.8};
  if(is.null(L)){L = as.numeric(quantile(x=z,probs=(1-p)/2))};
  if(is.null(U)){U = as.numeric(quantile(x=z,probs=(1+p)/2))};
  # Estimate parameters
  if(method=="CM"){theta = eNull.CM(z=z,f=f,L=L,U=U)};
  if(method=="ML"){theta = eNull.TN(z=z,L=L,U=U)};
  # Output
  return(theta);
}