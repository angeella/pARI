
# ARIpermutation
DOI: 10.5281/zenodo.3673779

Under construction.... 

<i class="fas fa-hammer"></i>


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
library(ARIpermutation)

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
alpha = 0.05
thr = 3.2

copes <- list()
sub_ids <- sapply(c(21:40),function(x) paste0(0,x))
for (sid in 1:length(sub_ids)) {  
  copes[[sid]] <- RNifti::readNifti(system.file("extdata", paste0("/sub-", sub_ids[sid] , ".nii.gz"), package = "ARIpermutation"))
  
}
mask <- system.file("extdata", "mask.nii.gz", package = "ARIpermutation")

out <- ARIpermCT(copes,thr=thr,mask=mask,alpha = alpha,family = "Simes")

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
