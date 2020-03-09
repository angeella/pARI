
# ARIpermutation
DOI: 10.5281/zenodo.3673779

**ARIpermutation** is the package developed to compute the All-Resolution Inference (ARI) method in the permutation framework. Therefore, this method doesn't 't assume any distribution about the null distribution of the p-values. It needs to satisfy the exchangeability assumption as all permutation-based methods.

As the parametric ARI, this method aims to compute simultaneous lower confidence bounds for the number of true discoveries, i.e., active voxels, in the fMRI framework. The function takes as input the list of copes, i.e., contrast maps, one for each subject, given by neuroimaging tools as FSL, SPM, etc. 

Having these data, you can insert a cluster map that can be the high-level output from some neuroimaging software, a region of interests (ROI), etc. If you want to construct these cluster maps using a supra-threshold statistic rule, you can specify the threshold into the argument ```thr``` of the function.

Therefore, the function ```ARIpermCT``` returns the lower bounds of true discoveries, i.e., active voxels, for each cluster coming from the cluster map inserted.

You can insert these cluster maps as many times as you want, because the Permutation-based ARI, as the parametric version, allows for **circular analysis**, still controlling for the multiplicity of inferences.

<!-- badges: start -->
<!-- badges: end -->

## Installation

You can install the released version of ARIpermutation with:

``` r
devtools::install_github("angeella/ARIpermutation")
```

## Simulation

Here, you can perform a toy example, using simulated data where the tests under the null hypotheses come from a Normal distribution with mean $0$ and variance $0.05$ and the tests under the alternative come from a Normal distribution with mean $10$ and variance $0.05$. We simulate $10$ tests under the null and $10$ under the alternative considering $10$ observations. Therefore, we have a matrix with dimensions $10 \times 20$, where the rows represent the observations and the columns the variables, i.e., tests.

We expect that the lower bound of the number of true discoveries, considering the full set of hypotheses, equals to $10$.

``` r
library(ARIpermutation)

m <- 20 #number of tests
n <- 10 $number of observations
X <- matrix(rnorm(0.5*m*n, 0, 0.05),ncol=n,nrow=0.5*m) #tests under the null
Y <- matrix(rnorm(0.5*m*n, 10, 0.05),ncol=n,nrow=0.5*m) #tests under the alternative
data <- cbind(X,Y) #full set of datasets

```
Then, we perform the sign-flipping test, using $2^n = 1024$ permutations, thanks to the function ```signTest```(type ```?ARIpermutation::signTest``` for more details):

``` r
pvalues <- signTest(data, 2^n)
```
and we plot it in $-log_{10}$ scale:

``` r
plot(-log(pvalues$pv,base = 10), pch = 20)
```
We create the p-values matrix where the rows represent the permutations and the columns the variables, i.e., the first row represents the observed p-values; the remain rows represent the p-values under the null distribution.

``` r
pv <- t(cbind(pvalues$pv,pvalues$pv_H0))
```

Then, we use the parametric approach considering the full set of hypotheses, i.e., ```ix``` equals ```c(1:10)```, using the function ```hommel``` and ```discoveries``` from the hommel package (type ```?hommel::hommel``` and ```hommel::discoveries``` for more details):

``` r
hom <- hommel(pv[1,], simes = TRUE)
discoveries(hom,ix = c(1:10),alpha = 0.1)

```
and the permutation-based one using the function ```SingleStepCT``` (type ```?ARIpermutation::SingleStepCT``` for more details)

``` r
SingleStepCT(data,ct = c(0,1),ix = c(1:10),alpha = 0.1,family = "Simes", B= 1000)[1]

```

We have at least $10$ true discoveries considering the full set of hypotheses.

## fMRI data 

This is a basic example using fMRI data from the [Auditory dataset](https://openneuro.org/datasets/ds000116/versions/00003). We need the list of copes:

``` r
copes <- list()
sub_ids <- sapply(c(21:40),function(x) paste0(0,x))
for (sid in 1:length(sub_ids)) {  
  copes[[sid]] <- RNifti::readNifti(system.file("extdata/AuditoryData", paste0("/sub-", sub_ids[sid] , ".nii.gz"), package = "ARIpermutation"))
  
}

```
the mask

``` r

mask <- system.file("extdata/AuditoryData", "mask.nii.gz", package = "ARIpermutation")

```
and the $\alpha$ level value and the threshold in order to perform the cluster map using a supra-threshold statistic rule: 


``` r
alpha = 0.1
thr = 3.2
```

then we can perform the Permutation-based ARI using the function ```ARIpermCT```(type ```?ARIpermutation::ARIpermCT``` for more details):

``` r
out <- ARIpermCT(copes,thr=thr,mask=mask,alpha = alpha,family = "Simes")
```

you can produce also the True Discovey Proportion brain map:

``` r
map_TDP(out,path= getwd(), name = "tdp", mask)
```

Then, you can compare it with the parametric method ARI using the [ARI](https://github.com/angeella/ARIbrain) package: 

``` r
Statmap <- system.file("extdata/AuditoryData", "Statmap.nii", package = "ARIpermutation")
mask <- system.file("extdata/AuditoryData", "mask.nii.gz", package = "ARIpermutation")
Pmap <- system.file("extdata/AuditoryData", "Pvaluemap.nii", package = "ARIpermutation")

#Create Clusters using a threshold equal to 3.2
Statmap = get_array(Statmap)
mask = get_array(mask)
Statmap[!mask]=0
clstr=cluster_threshold(Statmap>3.2)

res_ARI=ARIbrain::ARI(Pmap = Pmap, clusters= clstr, mask=mask, Statmap = Statmap)

```
# References

Rosenblatt, J. D., Finos, L., W., W. D., Solari, A., and Goeman, J. J. (2018). All-resolutions inference for brain imaging. NeuroImage, 181:786-796.

Hemerik, J., Solari, A., and Goeman, J. J. (2019). Permutation-based simultaneous confidence bounds for the false discovery proportion. Biometrika, 106(3):635-649.

Eklund, A., Nichols, E. T., and Knutsson, H. (2016). Cluster failure: Why fmri inferences for spatial extent have inflated false-positive rates. Pnas, 113(28):7900-7905.

Blanchard, G., Neuvial, P., and Roquain, E. (2019). Post hoc confidence bounds on false positives using reference families. Submitted to the Annals of Statistics.

# Did you find some bugs?

Please write to angela.andreella[\at]stat[\dot]unipd[\dot]it or insert a reproducible example using [reprex](https://github.com/tidyverse/reprex) on my [issue github page](https://github.com/angeella/ARIpermutation/issues).

