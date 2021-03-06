% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/signTest.R
\name{signTest}
\alias{signTest}
\title{Sign-Flipping Test}
\usage{
signTest(X, B = 1000, alternative = "two.sided", seed = NULL, mask = NULL, rand = F)
}
\arguments{
\item{X}{data where rows represents the variables and columns the observations}

\item{B}{by default \code{B = 1000}. Number of permutations.}

\item{alternative}{a character string referring to the alternative hypothesis, must be one of \code{"two.sided"} (default), \code{"greater"} or \code{"lower"}.}

\item{seed}{by default \code{seed=1234}. Integer value specifying the seed.}

\item{mask}{3D array of locicals (i.e. \code{TRUE/FALSE} in/out of the brain). Alternatively it may be a (character) NIfTI file name. If \code{mask=NULL}, it is assumed that non of the voxels have to be excluded.}

\item{rand}{by default \code{rand = FALSE}.}
}
\value{
Returns a list with the following objects: \code{Test} observed one sample t-test, \code{Test_H0} Test statistics under H0, \code{pv} observed p-values, \code{pv_H0} p-values under H0.
}
\description{
Performs sign-flipping, i.e. permutation, one-sample t-tests
}
\author{
Angela Andreella
}
