---
title: "README"
author: "Zachary McCaw"
date: "2018-09-25"
output: 
  html_document: 
    keep_md: TRUE
--- 

# Package Vignette




# Contents

* [Mixture Density](#mixture-density)
* [Null Density](#null-density)
* [Local FDR](#local-fdr)

# Mixture Density

## Density Estimation

Consider the following two-component normal mixture:

$$
\pi_{0}N(\mu_{0},\sigma_{0}^2) + (1-\pi_{0})N(\mu_{1},\sigma_{1}^{2})
$$
Below, $2\times 10^{3}$ observations are drawn from the two-component normal mixture with null proportion $\pi_{0}=0.95$, null parameters $(\mu_{0}=0,\sigma_{0}^{2}=1)$, and alternative parameters $(\mu_{1}=2,\sigma_{1}^{2}=1)$:


```r
library(LFDR);
set.seed(100);
# Data generation
D = rNM(n=2e3,k=2,pi=c(0.95,0.05),m=c(0,2),v=c(1,1));
# True component
Comp = D[,1];
# Observation
z = D[,2];
```

An exponential family density of order $d=5$ is fit to the sample by using 1. a natural spline basis, with knots on the sample quantiles, and 2. an orthogonal polynomial basis.  


```r
# Fit density via natural splines
fit.ns = Fit.EF(z=z,b=100,d=5,method="ns");
# Fit density via orthogonal polynomials
fit.poly = Fit.EF(z=z,b=100,d=5,method="poly");
```

<img src="README_files/figure-html/unnamed-chunk-4-1.png" style="display: block; margin: auto;" />

## Density Sampling

A sample $(z_{i}^{*})$ of size $n=10^{4}$ is drawn from fitted exponential family using the accept-reject (AR) method. The proposal distribution is normal, centered on the mean of the $(z_{i})$, with twice the standard deviation of the $(z_{i})$.


```r
# Sample fitted density
S = Sample.EF(z=z,f=fit.poly$f,m=mean(z),v=2*var(z),n=1e4,parallel=F);
```

The AR sample $(z_{i}^{*})$ allows for estimation of moments for the fitted exponential family. To match the generative model, the AR sample should have mean $2\pi_{0} = 0.1$ and variance $1+4\pi_{0}(1-\pi_{0}) \approx 1.2$:


```r
# Mean of AR sample
round(mean(S),digits=2);
# Variance of AR sample
round(var(S),digits=2);
```

# Empirical Null Distribution

## Estimation

The empirical null distribution is taken as normal, however the location and scale are not required to have the theoretical values of $\mu_{0}=0$ and $\sigma_{0}^{2}=1$. Rather, these parameters, and the proportion of observations $\pi_{0}$ arising from the null density, are estimated by central matching `CM` or maximum likelihood `ML`. 

Estimation requires specification of a *null neighborhood* $A_{0} = [L,U]$. Observations falling in the null neighborhood are assumed to have arisen from the null density. In central matching, $A_{0}$ is partitioned as $L = \xi_{0} < \cdots < \xi_{K} = U$. The log mixture density $\ln\hat{f}$ is evaluated over the partition points, and the evaluations $\ln\hat{f}_{i}$ are regressed on a quadratic function of the input. In maximum likelihood, $(\mu_{0},\sigma_{0}^{2})$ are estimated by fitting a truncated normal distribution to $A_{0}$. The null proportion is obtained from $(\mu_{0},\sigma_{0}^{2})$ and the proportion of $z$ scores falling in $A_{0}$. 


```r
# Estimation by central matching
M0 = eNull(z=z,L=-2,U=2,method="CM",f=fit.poly$f);
# Estimation by maximum likelihood
M1 = eNull(z=z,L=-2,U=2,method="ML");
# Tabulate estimates
nFit = round(rbind(as.numeric(M0),as.numeric(M1)),digits=2);
colnames(nFit) = c("Mean","Var","pi0");
print(nFit);
```

```
##      Mean  Var  pi0
## [1,] 0.06 1.12 0.98
## [2,] 0.05 1.07 0.97
```

# Local FDR

The local false discovery rate (lfdr) at $\zeta$ estimates the posterior probability that an observation with value $\zeta$ arose from the null component of the mixture: 

$$
\text{lfdr}(\zeta) = \frac{\pi_{0}f_{0}(\zeta)}{\pi_{0}f_{0}(\zeta) + (1-\pi_{0})f_{1}(\zeta)}
$$

The tail false discovery rate (tFDR) at $\zeta$ estimates the posterior probability that an observations with value as or more extreme than $\zeta$ arose from the null component of the mixture:

$$
\text{tFDR}(\zeta) = E[\text{fdr}(Z)|Z\geq \zeta] = \frac{\pi_{0}S_{0}(\zeta)}{\pi_{0}S_{0}(\zeta)+(1-\pi_{0})S_{1}(\zeta)}
$$

The function `FDR` estimates the local and tail false discovery rate functions using the observations `z`, a fitted mixture density $f(z)=\pi_{0}f_{0}(z)+(1-\pi_{0})f_{1}(z)$, and the parameters $(\mu_{0},\sigma_{0}^{2},\pi_{0})$ of the empirical null. Two functions are returned, `lfdr` and `tFDR` for calculating the local and tail false discovery rates, respectively. 


```r
# Estimate FDR curves
L = FDR(z=z,f=fit.poly$f,b=80,m=M1$m0,v=M1$v0,p0=M1$p0);
# Local FDR
lfdr = L$lfdr;
# Tail FDR
tFDR = L$tFDR;
# Evaluations for the observed data
l0 = lfdr(z);
t0 = tFDR(z);
```

<img src="README_files/figure-html/unnamed-chunk-9-1.png" style="display: block; margin: auto;" />
