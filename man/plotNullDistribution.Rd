% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotNullDistribution.R
\name{plotNullDistribution}
\alias{plotNullDistribution}
\title{Plot permutation p-values distribution}
\usage{
plotNullDistribution(P=NULL,family="simes",alpha = 0.05, 
path = getwd(), name = "plot", delta = 0,
copes=NULL,mask=NULL, alternative = "two.sided", rand = FALSE, B = 1000)
}
\arguments{
\item{P}{matrix of pvalues with dimensions \eqn{m \times B} used instead of the data matrix \code{X}. Default to @NULL.}

\item{family}{string character. Choose a family of confidence envelopes to compute the critical vector. 
from \code{"simes"}, \code{"aorc"}, \code{"beta"} and \code{"higher.criticism"}. If more than one critical vector are considered, it must be a vector. Default \code{"simes"}.}

\item{alpha}{numeric value in `[0,1]`. It expresses the alpha level to control the family-wise error rate. Default 0.1.}

\item{path}{character string. Path to save the plot. The path does not must end with \code{/}. Default to \code{getwd()}.}

\item{name}{character string. The name of file that will be used to save the plot. Default to name.}

\item{delta}{numeric value. It expresses the delta value, please see the references. 
Default to 0. If more than one critical vector are considered, \code{delta} must be a vector having length equals to the length of critical vectors specified in \code{family}.}

\item{copes}{list of NIfTI file if \code{P = NULL}. The list of copes, i.e., constrasts maps, one for each subject used to compute the statistical tests.}

\item{mask}{NIfTI file or character string. 3D array of logical values (i.e. \code{TRUE/FALSE} in/out of the brain). 
Alternatively it may be a (character) NIfTI file name. If \code{mask=NULL}, it is assumed that non of the voxels have to be excluded.}

\item{alternative}{character string. It refers to the alternative hypothesis, must be one of \code{"two.sided"} (default), \code{"greater"} or \code{"lower"}.}

\item{rand}{Boolean value. Default @FALSE. If \code{rand = TRUE}, the p-values are computed by \code{\link{rowRanks}}.}

\item{B}{numeric value. Number of permutations, default to 1000.}
}
\value{
Save a plot in \code{path} with name specified in \code{name} describing the p-values null distribution with critical value curve and observed pvalues in red.
}
\description{
Create a plot of permutation pvalues with corresponding specified critical vectors.
}
\examples{
\dontrun{
db <- simulateData(pi0 = 0.8, m = 100, n = 20, rho = 0)
out <- signTest(X = db)
pv <- cbind(out$pv, out$pv_H0)
plotNullDistribution(P = pv)
}

}
\author{
Angela Andreella
}
