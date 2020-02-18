
# ARIpermutation
DOI: 10.5281/zenodo.3673779

<!-- badges: start -->
<!-- badges: end -->

The aim of ARIpermutation is to ...

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
## Example

This is a basic example using fMRI data :

``` r

Statmap <- system.file("extdata", "Statmap.nii", package = "ARIpermutation")
mask <- system.file("extdata", "mask.nii.gz", package = "ARIpermutation")
Pmap <- system.file("extdata", "pv_par.nii", package = "ARIpermutation")
pvalues = system.file("extdata", "PvaluesPerm.Rda", package = "ARIpermutation")

alpha = 0.05
ct = c(0,1)
thr = 3.2
delta = 7

out1 <- ARIpermutation::ARIpermCT(Pmap, thr, mask=mask, Statmap= Statmap, alpha = alpha, pvalues = pvalues, ct = ct, family = "Simes",type="perm", delta = 7)

#Create Clusters using a threshold equal to 3.2
Statmap = get_array(Statmap)
mask = get_array(mask)
Statmap[!mask]=0
clstr=cluster_threshold(Statmap>3.2)

res_ARI=ARIbrain::ARI(Pmap = pvalue_name, clusters= clstr, mask=mask, Statmap = Statmap)

```
using the parametric method:

``` r

Statmap <- system.file("extdata", "Statmap.nii", package = "ARIpermutation")
mask <- system.file("extdata", "mask.nii.gz", package = "ARIpermutation")
Pmap <- system.file("extdata", "pv_par.nii", package = "ARIpermutation")
pvalues = system.file("extdata", "PvaluesPerm.Rda", package = "ARIpermutation")

#Create Clusters using a threshold equal to 3.2
Statmap = ARIbrain::get_array(Statmap)
mask = ARIbrain::get_array(mask)
Statmap[!mask]=0
clstr=cluster_threshold(Tmap>3.2)

res_ARI=ARIbrain::ARI(Pmap = pvalue_name, clusters= clstr, mask=mask, Statmap = Statmap)

```
