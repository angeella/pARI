
# ARIpermutation

<!-- badges: start -->
<!-- badges: end -->

The goal of ARIpermutation is to ...

## Installation

You can install the released version of ARIpermutation with:

``` r
devtools::install_github("angeella/ARIpermutation")
```

## Example

This is a basic example :

``` r
m <- 20
n <- 10
X <- matrix(runif(0.5*m*n, 0, 0.05),ncol=n,nrow=0.5*m)
Y <- matrix(runif(0.5*m*n, 1, 2),ncol=n,nrow=0.5*m)
data <- rbind(X,Y)
pvalues <- testByRandomization(data, 2^n)
plot(-log10(pvalues$p), pch = 20)
pv <- t(cbind(pvalues$p,pvalues$p0))
hom <- hommel(pv[1,], simes = TRUE)
discoveries(hom,ix = c(1:10),alpha = 0.1)
SingleStepCT(data,ct = c(0,1),ix = c(1:10),alpha = 0.1,family = "Simes", B= 1000)[1]

```

