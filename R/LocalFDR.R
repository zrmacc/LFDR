# Purpose: Calculate false discovery rates
# Updated: 180925

#' Local False Disocery Rate
#' 
#' Estimates the local false discovery rate from a collection of z scores. 
#' 
#' @param z Observations. 
#' @param f Mixture density. 
#' @param b Bins. 
#' @param m0 Mean of empirical null.
#' @param v0 Variance of empirical null.
#' @param p0 Null proportion. 
#' 
#' @return List including two functions, lfdr for calculating the local false
#'   discovery rate, and tFDR for calculating the tail false discovery rate.
#'
#' @importFrom stats approxfun integrate pnorm 
#' @export 

FDR = function(z,f=NULL,b=NULL,m0=0,v0=1,p0=1){
  # Bins
  if(is.null(b)){b=ceiling(length(z)/25)};
  # Obs
  n = length(z);
  # Range
  r = Support(z);
  R = as.numeric(r[2]-r[1]);
  # Bins
  bin = seq(from=r[1],to=r[2],length.out=(b+1));
  # Binwidth
  w = (R/b);
  # Midpoints
  mp = bin[1:(b)] + (1/2)*diff(bin);
  # Intervals
  A =  cbind(bin[1:b],bin[2:(b+1)]);
  # Discretize
  zd = cut(x=z,breaks=b,include.lowest=T,right=F);
  # Normal probability mass
  q = pnorm(bin,mean=m0,sd=sqrt(v0));
  q = q[2:(b+1)]-q[1:b];
  # Expected count
  e = n*(p0*q);
  # Smoothed, observed counts
  y = n*w*f(mp);
  # Local fdr
  aux = function(x){min(x,1)};
  lfdr = sapply(e/y,FUN=aux);
  # lfdr curve
  Curve = data.frame("zd"=levels(zd),"mp"=mp,"lfdr"=lfdr);
  # lfdr Step-function
  lfdr = approxfun(x=Curve$mp,y=Curve$lfdr,method="linear",f=0,rule=2);
  # Integration limits
  L = min(bin)-0.05*R;
  U = max(bin)+0.05*R;
  # Tail fdr probabilities
  i = NULL;
  tfdr = foreach(i=1:b,.combine=c) %do% {
    # Evaluation point
    x = mp[i];
    if(x<0){
      S0 = pnorm(q=x,mean=m0,sd=sqrt(v0));
      S = integrate(f=f,lower=L,upper=x)$value;
    } else {
      S0 = pnorm(q=x,mean=m0,sd=sqrt(v0),lower.tail=F);
      S = integrate(f=f,lower=x,upper=U)$value;
    }
    # Estimate tfdr
    Out = min(p0*S0/S,1);
    return(Out);
  }
  # Tail fdr function
  tFDR = approxfun(x=Curve$mp,y=tfdr,method="linear",f=0,rule=2);
  # Output
  Out = list("lfdr"=lfdr,"tFDR"=tFDR);
  return(Out);
}
