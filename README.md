
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
X <- matrix(rnorm(10*4000,5,sd=5),4000,10)

system.time(out1 <- testByRandomization(X,B=200))
```

