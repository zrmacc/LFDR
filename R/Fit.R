# Purpose: Fit an exponential family density
# Updated: 180923

#' Support
#' 
#' Return an interval containing the input data.
#' 
#' @param z Observations.
#' @param k Significant digits.
#' 
#' @return Numeric vector containing a lower and upper bound. 

Support = function(z,k=1){
  # Support
  m=min(z);
  M=max(z);
  # Coarsen
  s=10^(k);
  m = floor(m*s)/s;
  M = ceiling(M*s)/s;
  # Output
  Out = c(m,M);
  names(Out) = c("L","U");
  return(Out);
}

########################
# Estimation Function
########################

#' Fit Exponential Family Density
#' 
#' Fit an exponential family distribution of degree \code{d} to a 
#' sample of observations \code{z}. 
#' 
#' @param z Observations.
#' @param b Bins for discretizing
#' @param d Polynomial degree.
#' @param method Either "ns" for natural splines, or "poly" for orthogonal
#'   polynomials. Default is "poly".
#' 
#' @return A list containing a function \code{f} for evaluating the fitted
#'   density, and the estimated regression coefficients \code{beta}.
#' 
#' @importFrom splines ns
#' @importFrom stats coef glm.fit poisson poly quantile
#' @export 

Fit.EF = function(z,b=NULL,d=2,method="poly"){
  ## Input check
  if(is.null(b)){b=ceiling(length(z)/25)};
  # Method
  Choices = c("ns","poly");
  if(!(method%in%Choices)){stop("Select method from among: ns or poly.")};
  # Observed statistics
  n = length(z);
  # Support
  r = Support(z);
  # Bins
  bin = seq(from=r[1],to=r[2],length.out=(b+1));
  # Binwidth
  w = (r[2]-r[1])/b;
  # Midpoints
  mp = bin[1:(b)] + (1/2)*diff(bin);
  # Knot locations
  loc = quantile(z,probs=seq(from=0,to=d)/d);
  # Boundary knots
  loc.b = c(loc[1],loc[d+1]);
  # Interior knots
  loc.i = c(loc[2:d])
  # Discretize
  zd = cut(x=z,breaks=bin,include.lowest=T);
  y = as.numeric(table(zd));
  
  # Basis
  if(method=="ns"){
    # Spline matrix
    X = cbind(1,ns(mp,knots=loc.i,intercept=F,Boundary.knots=loc.b));
    colnames(X) = paste0(rep("s",times=(d+1)),seq(from=0,to=d));
  } else if(method=="poly"){
    # Polynomial matrix
    X = poly(mp,degree=d);
    # Extract parameters used to orthogonalize the polynomials
    a = attributes(X)$coefs;
    # Add intercept
    X = cbind(1,X);
    colnames(X) = paste0(rep("p",times=(d+1)),seq(from=0,to=d));
  };

  # Poisson regression
  off = rep(n*w,times=b);
  GLM = glm.fit(x=X,y=y,offset=log(off),family=poisson(link="log"));
  # Extract coefficients
  beta = coef(GLM);
  b0 = beta[1];
  b1 = beta[2:(d+1)];
  # Estimated density
  if(method=="ns"){
    f = function(x){
      p = splines::ns(x,knots=loc.i,intercept=F,Boundary.knots=loc.b);
      Out = as.numeric(exp(b0+MMP(p,b1)));
      return(Out);
    };
  } else if(method=="poly"){
    f = function(x){
      Out = as.numeric(exp(b0+MMP(poly(x,degree=d,coefs=a),b1)));
      return(Out);
    };
  };
  # Output
  Out = list("Method"=method,"f"=f,"beta"=beta);
  return(Out);
}

########################
# Sampling Function
########################

#' Sample Fitted Exponential Family
#' 
#' Sampling via the accept-reject algorithm. The proposal distribution
#' is taken as normal, with the option to specify the mean and variance. 
#' 
#' @param z Observations.
#' @param f Density to sample.
#' @param m Proposal mean.
#' @param v Proposal variance.
#' @param n Target number of sample to obtain.
#' @param parallel Run in parallel? Must register parallel backend first.
#' @param maxit Maximum number of accept-reject iterations.
#' 
#' @return Vector of realizations from the mixture density. 
#' 
#' @importFrom foreach '%do%' '%dopar%' foreach registerDoSEQ
#' @importFrom stats dnorm optimize rnorm runif var
#' @export

Sample.EF = function(z,f,m=NULL,v=NULL,n=1e3,parallel=F,maxit=25){
  # Support
  R = Support(z);
  # Proposal parameters
  if(is.null(m)){m=mean(z)};
  if(is.null(v)){v=var(z)};
  # Proposal density
  p = function(x){dnorm(x,mean=m,sd=sqrt(v))};
  # Likelihood ratio
  LR = function(x){f(x)/p(x)}
  # Search for optimum
  Opt = optimize(f=LR,interval=R,maximum=T);
  # Maximum likelihood ratio
  M = ceiling(Opt$objective*10)/10;
  if(!parallel){foreach::registerDoSEQ()};
  ## Implement accept reject
  Out = foreach(i=1:n,.combine=rbind) %dopar% {
    # Initialize parameters
    iter = 0;
    x = NA;
    keep = F;
    while(!keep){
      # Proposal
      z = rnorm(n=1,mean=m,sd=sqrt(v));
      # Uniform gague
      u = runif(n=1);
      # Accept-reject
      keep = (u < (1/M)*LR(z));
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