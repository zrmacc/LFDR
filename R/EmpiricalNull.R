# Purpose : Estimate empirical null

#' Central-Matching Estimate of the Empirical Null
#'
#' Estimate mean \eqn{\mu}, variance \eqn{\sigma^2}, and null proportion \eqn{\pi_{0}} 
#' by central matching.
#' 
#' @param z Observed z scores
#' @param f Mixture density
#' @param a Lower limit of null neighborhood
#' @param b Upper limit of null neighborhood
#' 
#' @importFrom stats coef lm.fit poly

eNullCm = function(z,f,a,b){
  # Evaluation grid
  x = seq(from=a,to=b,length.out=1e3);
  B = cbind(1,poly(x,degree=2,raw=T,simple=T));
  # Log mixture density
  y = log(f(x));
  # Regression log mixture density on polynomial in x
  lm0 = lm.fit(x=B,y=y);
  # Extract estimates
  beta0 = coef(lm0);
  # Variance
  v0 = as.numeric(-1/(2*beta0[3]));
  # Mean
  m0 = as.numeric(beta0[2]*v0);
  # Proportion null
  p0 = as.numeric(exp(beta0[1]+(1/2)*(m0^2/v0+log(2*pi*v0))));
  # Output
  Out = list("m0"=m0,"v0"=v0,"p0"=p0);
  return(Out);
}

#' Maximum Likelihood Estimation of the Empirical Null
#'
#' Estimate mean \eqn{\mu}, variance \eqn{\sigma^2}, and null proportion \eqn{\pi_{0}} 
#' by maximum likelihood.
#' 
#' @param z Observed z scores
#' @param a Lower limit of null neighborhood
#' @param b Upper limit of null neighborhood
#' @importFrom stats pnorm

eNullMl = function(z,a,b){
  # Total observations
  n = length(z);
  # Obs in null neighborhood
  z0 = z[z>=a & z<=b];
  n0 = length(z0);
  # Estimate theta
  theta = n0/n;
  # Fit truncated normal to obs in null neighborhood
  fit.tn = fitTruncNorm(z=z0,a=a,b=b,maxit=100); 
  # Location
  m0 = fit.tn$mu;
  # Scale
  s0 = fit.tn$sigma;
  # Null proportion
  H = function(m,s){pnorm((b-m)/s)-pnorm((a-m)/s)};
  p0 = theta/H(m0,s0);
  # Output
  Out = list("m0"=m0,"v0"=s0^2,"p0"=p0);
  return(Out);
}

#' Estimation of the Empirical Null
#' 
#' Estimate mean \eqn{\mu}, variance \eqn{\sigma^2}, and null proportion
#' \eqn{\pi_{0}}. Specify either the absolutely limits of the null neighborhood,
#' or the central proportion of the observations to serve as the null
#' neighborhood.
#' 
#' @param z Observed z scores
#' @param a Lower limit of null neighborhood
#' @param b Upper limit of null neighborhood
#' @param p Central proportion of z scores belonging to null neighborhood
#' @param method Either 1. CM central matching, or 2. ML maximum likelihood
#' @param f Estimated mixture density if \code{method="CM"}.
#' 
#' @importFrom stats quantile
#' @export 

eNull = function(z,a,b,p,method="ML",f){
  # Check inputs
  if(!(method %in% c("CM","ML"))){stop("Select method from CM or ML.")};
  if(method=="CM"&missing(f)){stop("For central matching, provide an estimate of the mixture density.")};
  if(method=="CM"&missing(p)){p=0.5};
  if(method=="ML"&missing(p)){p=0.8};
  if(missing(a)){a = quantile(x=z,probs=(1-p)/2);};
  if(missing(b)){b = quantile(x=z,probs=(1+p)/2);};
  # Estimate parameters
  if(method=="CM"){theta = eNullCm(z=z,f=f,a=a,b=b);};
  if(method=="ML"){theta = eNullMl(z=z,a=a,b=b)};
  # Output
  return(theta);
}