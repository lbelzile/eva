% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gpdImPb.R
\name{gpdImPb}
\alias{gpdImPb}
\title{GPD Bootstrapped Information Matrix (IM) Test}
\usage{
gpdImPb(data, inner, outer, allowParallel = FALSE, numCores = 1)
}
\arguments{
\item{data}{Data should be in vector form.}

\item{inner}{Number of bootstrap replicates for the covariance estimate.}

\item{outer}{Number of bootstrap replicates for critical values.}

\item{allowParallel}{Should the outer bootstrap procedure be run in parallel or not. Defaults to false.}

\item{numCores}{If allowParallel is true, specify the number of cores to use.}
}
\value{
\item{statistic}{Test statistic.}
\item{p.value}{P-value for the test.}
\item{theta}{Estimate of theta for the initial dataset.}
\item{effective_bootnum}{Effective number of outer bootstrap replicates used (only those that converged are used).}
}
\description{
Runs the IM Test using a two-step iterative procedure, to boostrap the covariance estimate and critical values. See reference for details.
}
\details{
Warning: This test can be very slow, since the covariance estimation is nested within the outer replicates. It would be
recommended to use a small number of replicates for the covariance estimate (at most 50).
}
\examples{
## Not run
# x <- rgpd(200, loc = 0, scale = 1, shape = 0.2)
# gpdImPb(x, inner = 20, outer = 99)
}
\references{
Dhaene, G., & Hoorelbeke, D. (2004). The information matrix test with bootstrap-based covariance matrix estimation. Economics Letters, 82(3), 341-347.
}
