# Purpose : Generate Z statistics from a mixture of normal densities

#' Generate Z's from Two-Component Normal-Mixutre
#' 
#' @param pi0 Proportion of Z's from the null density.
#' @param mu1 Location of alternative density.
#' @param sigma1 Spread of alternative density.
#' @param n Number of statistics to generate.
#' @param symm If TRUE, the alternative is symmetric about zero.
#' 
#' @importFrom stats rbinom rnorm 
#' @export 

genMixZ = function(pi0,mu1,sigma1=1,n,symm=F){
  # Number of null z
  n0 = round(pi0*n);
  # Number of alternative z
  n1 = n-n0;
  # Null z
  z.n = rnorm(n=n0);
  # Alternative z
  z.a = rnorm(n=n1,mean=mu1,sd=sigma1);
  # Symmetrize alternative
  if(symm){
    s.a = rbinom(n=n1,size=1,prob=0.5);
    s.a = 2*(s.a)-1;
    z.a = z.a*(s.a);
  }
  # Output
  z = c(z.n,z.a);
  Out = cbind("z"=z,"Density"=c(rep(0,n0),rep(1,n1)));
  # Data
  return(Out);
}