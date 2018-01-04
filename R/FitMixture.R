# Purpose : Fit an exponential family density to Z scores 

#' Range
#' 
#' Report an interval containing the input data.
#' @param z Data.
#' @param k Significant digits.
#' @export 

Range = function(z,k=1){
  # Support
  m=min(z);
  M=max(z);
  s=exp(k*log(10));
  # Coarsen
  m = floor(m*s)/s;
  M = ceiling(M*s)/s;
  return(c(m,M));
}

#' Fit Natural Spline Basis Exponential Family 
#' 
#' @param z Normalized statistics.
#' @param k Bins.
#' @param l Polynomial degree.
#' 
#' @importFrom splines ns
#' @importFrom stats coef glm.fit poisson quantile
#' @export

fitNsDens = function(z,k,l){
  # Observed statistics
  n = length(z);
  # Range
  r = Range(z);
  # Bins
  b = seq(from=r[1],to=r[2],length.out=(k+1));
  # Binwidth
  d = (r[2]-r[1])/k;
  # Midpoints
  mp = b[1:(k)] + (1/2)*diff(b);
  # Knot locations
  loc = quantile(z,probs=seq(from=0,to=l)/l);
  loc.b = c(loc[1],loc[l+1]);
  loc.i = c(loc[2:l])
  # Discretize
  zd = cut(x=z,breaks=b,include.lowest=T)
  # Design matrix
  A = splines::ns(mp,knots=loc.i,intercept=F,Boundary.knots=loc.b);
  C = cbind("y"=as.numeric(table(zd)),1,A);
  colnames(C)[2:(2+l)] = paste0(rep("s",times=l),seq(from=0,to=l));
  # Poisson regression
  off = rep(n*d,times=k);
  GLM = glm.fit(x=C[,2:ncol(C)],y=C[,1],offset=log(off),family=poisson(link="log"));
  beta = coef(GLM);
  b0 = beta[1];
  b1 = beta[2:(l+1)];
  # Return estimated density
  f = function(x){
    p = splines::ns(x,knots=loc.i,Boundary.knots=loc.b);
    return(as.numeric(exp(b0+p %*% b1)));
  };
  f = Vectorize(f);
  # Output
  Out = list("f"=f,"beta"=beta);
  return(Out);
}

#' Fit Polynomial Basis Exponential Family 
#' 
#' @param z Normalized statistics.
#' @param k Bins.
#' @param l Polynomial degree.
#' 
#' @importFrom stats coef glm.fit poisson poly quantile
#' @export

fitPolyDens = function(z,k,l){
  # Observed statistics
  n = length(z);
  # Range
  r = Range(z);
  # Bins
  b = seq(from=r[1],to=r[2],length.out=(k+1));
  # Binwidth
  d = (r[2]-r[1])/k;
  # Midpoints
  mp = b[1:(k)] + (1/2)*diff(b);
  # Discretize
  zd = cut(x=z,breaks=b,include.lowest=T)
  # Design matrix
  A = poly(mp,degree=l);
  a = attributes(A)$coefs; # Extracts parameters used to orthogonalize the polynomials
  C = cbind("y"=as.numeric(table(zd)),1,A);
  colnames(C)[2:(2+l)] = paste0(rep("x",times=l),seq(from=0,to=l));
  # Poisson regression
  off = rep(n*d,times=k);
  GLM = glm.fit(x=C[,2:ncol(C)],y=C[,1],offset=log(off),family=poisson(link="log"));
  beta = coef(GLM);
  b0 = beta[1];
  b1 = beta[2:(l+1)]
  # Return estimated density
  f = function(x){
    return(as.numeric(exp(b0+poly(x,degree=l,coefs=a) %*% b1)));
  };
  f = Vectorize(f);
  # Output
  Out = list("f"=f,"beta"=beta);
  return(Out);
}

#' Sample Fitted Exponential Family
#' 
#' Sampling via the accept-reject algorithm. The proposal distribution
#' is taken as normal, with the option to specify the mean and variance. 
#' 
#' @param f Density to sample. 
#' @param m Proposal mean.
#' @param s Proposal standard deviation.
#' @param n Approximate number of sample to obtain. 
#' @param R Range of data.
#' @param parallel Run in parallel? 
#' @param maxit Maximum iterations of accept-reject.
#' 
#' @importFrom foreach '%do%' '%dopar%' foreach registerDoSEQ
#' @importFrom stats dnorm optimize rnorm runif
#' @export

sampleFit = function(f,m=0,s=1,n,R,parallel=F,maxit=25){
  # Proposal density
  p = function(x){dnorm(x,mean=m,sd=s)};
  # Likelihood ratio
  LR = function(x){f(x)/p(x)}
  # Log likelihood ratio
  lLR = function(x){log(f(x))-log(p(x))};
  # Search for optimum
  Opt = optimize(f=LR,interval=R,maximum=T);
  # Maximum likelihood ratio
  k = ceiling(Opt$objective*10)/10;
  if(!parallel){foreach::registerDoSEQ()};
  ## Implement accept reject
  Out = foreach(i=1:n,.combine=rbind) %dopar% {
    # Initialize parameters
    iter = 0;
    x = NA;
    keep = F;
    while(!keep){
      # Proposal
      z = rnorm(n=1,mean=m,sd=s);
      # Uniform gague
      u = runif(n=1);
      # Accept-reject
      keep = (u < (1/k)*LR(z));
      # If accepted, break
      if(keep){x = z; break};
      # If rejected, increment
      iter = iter + 1;
      if(iter>=maxit){break};
    } # End while
    return(x)
  }
  # Output
  Out = as.numeric(Out[!is.na(Out)]);
  return(Out);
}