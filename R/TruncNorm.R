#' @useDynLib LFDR, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' Initial Estimates of Mean and Scale for Truncated Normal
#' 
#' @param z Data
#' @param a Lower truncation limit
#' @param b Upper truncation limit
#' @param maxit Maximum iterations
#' 
#' @importFrom stats median pnorm quantile qnorm uniroot var

initEst = function(z,a,b,maxit=5){
  # Auxiliary functions
  alpha = function(m,s){(a-m)/s};
  beta = function(m,s){(b-m)/s};
  H = function(m,s){pnorm(beta(m,s))-pnorm(alpha(m,s))};
  # Initialize scale
  q = quantile(z,probs=c(0.025,0.975));
  s0 = as.numeric((q[2]-q[1])/2);
  # Sample moments
  M = median(z);
  V = var(z);
  r = q[2]-q[1];
  # Iterate between estimation of m and s
  for(i in 1:maxit){
    # Method of moments for m
    g = function(m){
      Z = H(m,s0);
      e1 = (pnorm(alpha(m,s0)) + pnorm(beta(m,s0)))/2;
      E = m + s0*qnorm(e1)
      return(M-E);
    }
    m0 = uniroot(f=g,lower=a,upper=b)$root;
    # Method of moments for s
    h = function(s){
      Z = H(m0,s);
      e1 = (alpha(m0,s)*dnorm(alpha(m0,s))-beta(m0,s)*dnorm(beta(m0,s)))/Z;
      e2 = ((dnorm(alpha(m0,s))-dnorm(beta(m0,s)))/Z)^2;
      E = (s^2)*(1+e1-e2);
      return(V-E);
    }
    s0 = uniroot(f=h,lower=0.1,upper=r)$root;
  }
  # Output
  theta0 = c("m0"=m0,"s0"=s0);
  return(theta0);
}

#' Estimate Meand and Scale of Truncated Normal
#' 
#' @param z Data.
#' @param a Lower truncation limit.
#' @param b Upper truncation limit.
#' @param maxit Maximum iterations.
#' @param trace Trace log likelihood.
#' 
#' @importFrom stats dnorm pnorm
#' @export 

fitTruncNorm = function(z,a,b,maxit=25,trace=F){
  # Observations
  n = length(z);
  # Auxiliary functions
  phi = dnorm;
  alpha = function(m,s){(a-m)/s};
  beta = function(m,s){(b-m)/s};
  # Renormalizing constant, and derivatives
  H = function(m,s){pnorm(beta(m,s))-pnorm(alpha(m,s))};
  Hm = function(m,s){-(1/s)*(phi(beta(m,s))-phi(alpha(m,s)))}
  Hs = function(m,s){-(1/s)*(beta(m,s)*phi(beta(m,s))-alpha(m,s)*phi(alpha(m,s)))}
  Hmm = function(m,s){-(1/s)^2*(beta(m,s)*phi(beta(m,s))-alpha(m,s)*phi(alpha(m,s)))}
  Hms = function(m,s){-(1/s)^2*(beta(m,s)^2*phi(beta(m,s))-alpha(m,s)^2*phi(alpha(m,s)))}
  Hss = function(m,s){-(1/s)^2*(beta(m,s)^3*phi(beta(m,s))-alpha(m,s)^3*phi(alpha(m,s)))}
  # Score equations
  Um = function(m,s){-(n/H(m,s))*Hm(m,s)+(1/s)^2*sum(z-m)};
  Us = function(m,s){-(n/s)-(n/H(m,s))*Hs(m,s)+(1/s)^3*sum((z-m)^2)};
  # Reparameterize the scale to ensure positivity
  Ueta = function(m,eta){s = exp(eta);return(Us(m,s)*s)};
  U = function(theta){
    m = theta[1];
    eta = theta[2];
    s = exp(eta);
    return(c(Um(m,s),Ueta(m,eta)));
  }
  # Information components
  Imm = function(m,s){-(n/H(m,s)^2)*(Hm(m,s))^2+(n/H(m,s))*Hmm(m,s)+n/(s^2)}
  Ims = function(m,s){-(n/H(m,s)^2)*(Hm(m,s)*Hs(m,s))+(n/H(m,s))*Hms(m,s)}
  Imeta = function(m,eta){
    s = exp(eta);
    return(Ims(m,s)*s);
  }
  Iss = function(m,s){-(n/H(m,s)^2)*(Hs(m,s))^2+(n/H(m,s))*Hss(m,s)+2*n/s^2}
  Ietaeta = function(m,eta){
    s = exp(eta);
    return(Iss(m,s)*s^2);
  }
  Info = function(theta){
    m = theta[1]; 
    eta = theta[2];
    s = exp(eta);
    return(matrix(c(Imm(m,s),Imeta(m,eta),Imeta(m,eta),Ietaeta(m,eta)),byrow=T,nrow=2));
    };
  # Objective function
  Q = function(theta){
    m = theta[1];
    s = exp(theta[2]);
    return(-n*log(s)-n*log(H(m,s))-1/(2*s^2)*sum((z-m)^2));
  }
  # Method of moments estimators
  theta0 = initEst(z,a,b,maxit=5);
  theta0[2] = log(theta0[2]);
  # Initial objective
  q00 = q0 = Q(theta0);
  # Initial likelihood
  if(trace){
    cat(paste0("Initial Likelihood: ",round(q0,digits=2)));
    cat("\n");
    };
  # Iterate
  for(i in 1:maxit){
    # Score
    s0 = U(theta0);
    # Information
    i0 = Info(theta0);
    # NR Update
    d = Delta(Info=i0,S=s0);
    # NR update
    theta1 = theta0 + d;
    # New objective
    q1 = Q(theta1);
    # Indicator objective decreased
    if(q1<=q0){break};
    theta0 = theta1;
    q0 = q1;
  }
  # Final likelihood
  q2 = q0;
  # Update
  if(trace){
    cat(paste0("Final Likelihood: ",round(q2,digits=2)));
    cat("\n");
    cat(paste0("Overall Likelihood Improvement: ",signif(q2-q00,digits=2)));
    cat("\n");
    };
  # Final estimates
  names(theta0) = NULL;
  Out = list("mu"=theta0[1],"sigma"=exp(theta0[2]));
  return(Out);
}