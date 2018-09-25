# Purpose: Estiated parameters for truncated normal distribution
# Updated: 180924

########################
# Renormalizing Constant
########################

# Verified numerically. 

#' Truncated Normal Renormalizing Constant
#' 
#' Calculates the renormalizing constant of the truncated
#' normal distribution, and its derivatives.
#' 
#' @param m Location parameter.
#' @param s Scale parameter.
#' @param L Lower truncation limit.
#' @param U Upper truncation limit.
#' @param k Order.
#' @param dirn Direction.
#' 
#' @importFrom stats dnorm pnorm

H = function(m,s,L,U,k=0,dirn){
  # Alpha and beta
  a = (L-m)/s;
  b = (U-m)/s;
  # Initialize
  Out = 0;
  # Order
  if(k==0){
  # Zeroth derivative
    Out = pnorm(b)-pnorm(a);
  } else if(k==1){
  ## First derivatives 
    # In location
    if(dirn=="m"){
      Out = -(1/s)*(dnorm(b)-dnorm(a));
    };
    # In scale
    if(dirn=="s"){
      Out = -(1/s)*(b*dnorm(b)-a*dnorm(a));
    };
  } else if(k==2){
  ## Second derivatives
    # In location
    if(dirn=="m"){
      Out = -(1/s)^2*(b*dnorm(b)-a*dnorm(a));
    };
    # In scale
    if(dirn=="s"){
      Out = -(1/s)^2*(b*(b^2-2)*dnorm(b)-a*(a^2-2)*dnorm(a));
    }
    # Mixed partial
    if(dirn=="ms"){
      Out = -(1/s)^2*((b^2-1)*dnorm(b)-(a^2-1)*dnorm(a));
    }
  };
  # Output
  return(Out);
};

########################
# Log Likelihood
########################

#' Truncated Normal Log Likelihood
#' 
#' Evaluates the sample log likelihood of the truncated
#' normal distribution. 
#' 
#' @param z Observations.
#' @param m Location parameter.
#' @param t Log of scale parameter.
#' @param L Lower truncation limit.
#' @param U Upper truncation limit.

qTN = function(z,m,t,L,U){
  # Scale parameter
  s = exp(t);
  # Obs
  n = length(z);
  # Normalizing constant
  H0 = H(m=m,s=s,L=L,U=U,k=0);
  # Objective
  Out = -n*log(s)-n*log(H0)-sum((z-m)^2)/(2*s^2);
  # Output
  return(Out);
}

########################
# Score Equations
########################

# Verified numerically. 

#' Truncated Normal Score Equations
#' 
#' Evaluates the gradient of the sample log likelihood. 
#' 
#' @param z Observations.
#' @param m Location parameter.
#' @param t Log of scale parameter.
#' @param L Lower truncation limit.
#' @param U Upper truncation limit.

uTN = function(z,m,t,L,U){
  # Obs
  n = length(z);
  # Scale parameter
  s = exp(t);
  # Normalizing constant
  H0 = H(m=m,s=s,L=L,U=U,k=0);
  # Derivative in m
  Hm = H(m=m,s=s,L=L,U=U,k=1,dirn="m");
  # Derivative in s
  Hs = H(m=m,s=s,L=L,U=U,k=1,dirn="s");
  # Score for location
  um = -(n/H0)*Hm+sum((z-m))/(s^2);
  # Score for log scale
  ut = (-(n/s)-(n/H0)*Hs+sum((z-m)^2)/(s^3))*s;
  # Format
  Out = c(um,ut);
  names(Out) = c("um","ut");
  # Output
  return(Out);
}

########################
# Information Matrix
########################

# Verified numerically. 

#' Truncated Normal Observed Information
#' 
#' Evaluates the observed information of the sample log likelihood. 
#' 
#' @param z Observations.
#' @param m Location parameter.
#' @param t Log of scale parameter.
#' @param L Lower truncation limit.
#' @param U Upper truncation limit.

infoTN = function(z,m,t,L,U){
  # Obs
  n = length(z);
  # Scale parameter
  s = exp(t);
  # Normalizing constant
  H0 = H(m=m,s=s,L=L,U=U,k=0);
  # Derivatives in m
  Hm  = H(m=m,s=s,L=L,U=U,k=1,dirn="m");
  Hmm = H(m=m,s=s,L=L,U=U,k=2,dirn="m");
  # Derivatives in s
  Hs  = H(m=m,s=s,L=L,U=U,k=1,dirn="s");
  Hss = H(m=m,s=s,L=L,U=U,k=2,dirn="s");
  # Mixed partial
  Hms = H(m=m,s=s,L=L,U=U,k=2,dirn="ms");
  # Information for m
  Imm = -n/(H0^2)*(Hm^2)+(n/H0)*Hmm+n/(s^2);
  # Cross information
  Ims = -n/(H0^2)*(Hm*Hs)+n/H0*Hms+2*sum((z-m))/(s^3);
  Imt = Ims*s;
  # Gradient
  ut = uTN(z=z,m=m,t=t,L=L,U=U)[2];
  # Information for t
  hss = n/(s^2)+n/(H0^2)*(Hs^2)-n/H0*(Hss)-3*sum((z-m)^2)/(s^4);
  Itt = -1*(hss*s^2+ut);
  # Output
  Out = matrix(c(Imm,Imt,Imt,Itt),nrow=2);
  rownames(Out) = colnames(Out) = c("m","t");
  return(Out);
}

########################
# Fit Truncated Normal
########################

#' Estimate Mean and Scale of Truncated Normal
#' 
#' Estimates the location and scale of a truncated normal distribution. 
#' 
#' @param z Observations.
#' @param L Lower truncation limit.
#' @param U Upper truncation limit.
#' @param maxit Maximum iterations.
#' @param eps Tolerance.
#' @param report Report fitting progress?
#' 
#' @return Vector containing the estimated location and scale of the truncated
#'   normal distribution.
#' 
#' @importFrom stats dnorm pnorm
#' @export 

fitTN = function(z,L,U,maxit=25,eps=1e-8,report=F){
  
  # Update function
  Update = function(theta){
    # Current parameters
    m0 = theta$m;
    t0 = theta$t;
    # Initial log likelihood
    q0 = qTN(z=z,m=m0,t=t0,L=L,U=U);
    # Score
    u0 = uTN(z=z,m=m0,t=t0,L=L,U=U);
    # Information
    J0 = infoTN(z=z,m=m0,t=t0,L=L,U=U);
    # Proposal
    prop = c(m0,t0)+as.numeric(MMP(matInv(J0),u0));
    m1 = prop[1];
    t1 = prop[2];
    # Final log likelihood
    q1 = qTN(z=z,m=m1,t=t1,L=L,U=U);
    # Increments
    d = q1-q0;
    # Output
    Out = list("m"=m1,"t"=t1,"d"=d);
  }
  
  # Initialize at zero
  theta0 = list("m"=0,"t"=0);
  ## Maximzation
  for(i in 1:maxit){
    # Update
    theta1 = Update(theta0);
    # Accept if increment is positive
    if(theta1$d>0){
      theta0 = theta1;
      if(report){cat("Objective increment: ",signif(theta1$d,digits=3),"\n")}
    }
    # Terminate if increment is below tolerance
    if(theta1$d<eps){
      rm(theta1);
      break;
    }
  };
  
  # Final estimates
  Out = c(theta0$m,exp(theta0$t));
  names(Out) = c("m","s");
  return(Out);
}
