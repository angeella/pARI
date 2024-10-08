% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dI.R
\name{dI}
\alias{dI}
\title{Lower bound for the number of true discoveries}
\usage{
dI(ix, cv, pvalues, iterative, approx, ncomb, ...)
}
\arguments{
\item{ix}{Numeric vector: set-wise hypotheses considered.}

\item{cv}{Numeric vector: critical vector computed by \code{\link{criticalVector}}.}

\item{pvalues}{If \code{iterative = TRUE} you must put here the matrix of \eqn{p}-values with 
dimensions \eqn{m \times B} where \eqn{m} is the number of variables and \eqn{B} the number of permutations. 
Instead, if \code{iterative = FALSE}, you can put directly the vector of \eqn{m} observed \eqn{p}-values.}

\item{iterative}{Boolean value. If \code{iterative = TRUE}, the iterative method is applied (computationally demanding). Default to \code{FALSE}. Please see the reference below.}

\item{approx}{Boolean value. Default to \code{TRUE}. If you are analyzing high dimensional data, we suggest to put \code{approx = TRUE} to speed up the computation time. Please see the reference below.}

\item{ncomb}{Numeric value. If \code{approx = TRUE}, you must decide how many random sub collections (level of approximation) considered. Default to 100.}

\item{...}{Further arguments for the iterative approach, i.e., \code{iterative = TRUE}.}
}
\value{
Numeric value: the lower confidence bound for the number of true discoveries concerning the cluster \code{ix} specified.
}
\description{
Calculates 1-\eqn{\alpha} lower confidence bound for the set-wise of false null hypotheses.
}
\examples{
db <- simulateData(pi0 = 0.7, m = 100, n = 20, rho = 0)
out <- signTest(X = db)
pv <- cbind(out$pv, out$pv_H0)
cv <- criticalVector(pvalues = pv, family = "simes", lambda = 0.1, alpha = 0.1)
dI(ix = c(1:100), cv = cv, pvalues = pv)
}
\references{
Andreella, A., Hemerik, J., Finos, L., Weeda, W., & Goeman, J. (2023). Permutation-based true discovery proportions for functional magnetic resonance imaging cluster analysis. Statistics in Medicine, 42(14), 2311-2340.
}
\author{
Angela Andreella
}
