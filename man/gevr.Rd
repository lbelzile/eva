% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gevr.R
\name{gevr}
\alias{gevr}
\alias{dgevr}
\alias{rgevr}
\alias{qgev}
\alias{pgev}
\title{The GEVr Distribution}
\usage{
dgevr(x, loc = 0, scale = 1, shape = 0, log.d = FALSE)

rgevr(n, r, loc = 0, scale = 1, shape = 0)

qgev(p, loc = 0, scale = 1, shape = 0, lower.tail = TRUE, log.p = FALSE)

pgev(q, loc = 0, scale = 1, shape = 0, lower.tail = TRUE, log.p = FALSE)
}
\arguments{
\item{x}{Vector or matrix of observations. If x is a matrix, each row is taken to be a new observation.}

\item{loc, scale, shape}{Location, scale, and shape parameters. Can be vectors, but
the lengths must be appropriate.}

\item{log.d}{Logical: Whether or not to return the log density. (FALSE by default)}

\item{n}{Number of observations}

\item{r}{Number of order statistics for each observation.}

\item{p}{Vector of probabilities.}

\item{lower.tail}{Logical: If TRUE (default), probabilities are P[X <= x] otherwise, P[X > x].}

\item{log.p}{Logical: If TRUE, probabilities p are given as log(p). (FALSE by default)}

\item{q}{Vector of quantiles.}
}
\description{
Random number generation (rgevr) and density (dgevr) functions for the GEVr distribution with parameters loc, scale, and shape.
Also, quantile function (qgev) and cumulative distribution function (pgev) for the GEV1 distribution.
}
\details{
GEVr data (in matrix x) should be of the form \eqn{x[i,1] > x[i, 2] > \cdots > x[i, r]} for each observation
\eqn{i = 1, \ldots, n}. Note that currently the quantile and cdf functions are only for the GEV1 distribution. The GEVr
distribution is also known as the r-largest order statistics model and is a generalization of the block maxima model (GEV1).
The density function is given by \deqn{f_r (x_1, x_2, ..., x_r | \mu, \sigma, \xi) = \sigma^{-r} \exp\Big\{-(1+\xi z_r)^{-\frac{1}{\xi}}
- \left(\frac{1}{\xi}+1\right)\sum_{j=1}^{r}\log(1+\xi z_j)\Big\}} for some location parameter \eqn{\mu},
scale parameter \eqn{\sigma > 0}, and shape parameter \eqn{\xi}, where \eqn{x_1 > \cdots > x_r}, \eqn{z_j = (x_j - \mu) / \sigma},
and \eqn{1 + \xi z_j > 0} for \eqn{j=1, \ldots, r}. When \eqn{r = 1}, this distribution is exactly the GEV distribution.
}
\examples{
## Plot the densities of the heavy and bounded upper tail forms of GEVr
set.seed(7)
dat1 <- rgevr(1000, 1, loc = 0, scale = 1, shape = 0.25)
dat2 <- rgevr(1000, 1, loc = 0, scale = 1, shape = -0.25)
hist(dat1, col = rgb(1, 0, 0, 0.5), xlim = c(-5, 10), ylim = c(0, 0.4),
     main = "Histogram of GEVr Densities", xlab = "Value", freq = FALSE)
hist(dat2, col = rgb(0, 0,1, 0.5), add = TRUE, freq = FALSE)
box()

## Generate sample with decreasing trend in location parameter
x <- rgevr(10, 2, loc = 10:1, scale = 1, shape = 0.1)
dgevr(x, loc = 10:1, scale = 10:1, shape = 0.1)

## Incorrect parameter specifications
# rgevr(10, 2, loc = 5:8, scale = 1, shape = 0.1)
# rgevr(1, 2, loc = 5:8, scale = 1:2, shape = 0.1)
}
\references{
Coles, S. (2001). An introduction to statistical modeling of extreme values (Vol. 208). London: Springer.
}
