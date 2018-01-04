# Purpose : Calculate local FDR

#' Local False Disocery Rate
#' 
#' Calculate the local false discovery rate for a set of z scores. In addition, returns 
#' funcations for calculating the local and tail are FDRs.
#' 
#' @param z Observed z scores
#' @param f Mixture density function
#' @param k Number of bins
#' @param m Mean of empirical null
#' @param v Variance of empirical null
#' @param p0 Null proportion
#' 
#' @return List including 1. a data frame, 2. a function for calculating the
#'   local FDR, and 3. a function for calculating the tail area FDR. The data
#'   frame includes the observed z scores, the discretized interval into which
#'   the z score fell, the midpoint of the interval, and the estimated local FDR
#'   for that interval.
#'
#' @importFrom stats approxfun integrate pnorm 
#' @export 

LFDR = function(z,f,k,m=0,v=1,p0=1){
  # Obs
  n = length(z);
  # ID to recover original order of observations
  ID = seq(1:n);
  # Range
  r = Range(z);
  # Bins
  b = seq(from=r[1],to=r[2],length.out=(k+1));
  # Binwidth
  d = (r[2]-r[1])/k;
  # Midpoints
  mp = b[1:(k)] + (1/2)*diff(b);
  # Intervals
  A =  cbind(b[1:k],b[2:(k+1)]);
  # Discretize
  zd = cut(x=z,breaks=b,include.lowest=T,right=F);
  Map = data.frame(ID,z,zd);
  colnames(Map) = c("ID","z","zd");
  # Normal probability mass
  q = pnorm(b,mean=m,sd=sqrt(v));
  q = q[2:(k+1)]-q[1:k];
  # Expected count
  e = n*(p0*q);
  # Smoothed, observed counts
  y = n*d*f(mp);
  # Local fdr
  aux = function(x){min(x,1)};
  lfdr = sapply(e/y,FUN=aux);
  # lfdr curve
  Curve = data.frame("zd"=levels(zd),"mp"=mp,"lfdr"=lfdr);
  # lfdr Step-function
  lfdr = approxfun(x=b,y=c(Curve$lfdr,Curve$lfdr[nrow(Curve)]),method="constant",f=0,rule=2);
  # lfdr mappings
  Map = merge(x=Map,all.x=T,y=Curve,by="zd");
  Map = Map[,c("ID","z","zd","mp","lfdr")];
  # Remove original order
  Map = Map[order(Map$ID),];
  Map = Map[,2:5]
  rownames(Map) = NULL;
  # Tail fdr function
  tFDR = function(x){
    g = function(x){f(x)*lfdr(x)};
    if(x<=0) {
      a = integrate(g,lower=r[1],upper=x,abs.tol=1e-4,stop.on.error=F)$value;
      b = integrate(f,lower=r[1],upper=x,abs.tol=1e-4,stop.on.error=F)$value;
    } else {
      a = integrate(g,lower=x,upper=r[2],abs.tol=1e-4,stop.on.error=F)$value;
      b = integrate(f,lower=x,upper=r[2],abs.tol=1e-4,stop.on.error=F)$value;
    };
    return(a/b);
  }
  tFDR = Vectorize(tFDR);
  # Output
  Out = list("Map"=Map,"lfdr"=lfdr,"tFDR"=tFDR);
  return(Out);
}