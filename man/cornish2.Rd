% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cornish.R
\name{cornish2}
\alias{cornish2}
\title{Calculate the Cornish-Fisher two-dimensional quantile approximation}
\usage{
cornish2(alpha, klist, nn = 1, xv = NULL, prestd = FALSE)
}
\arguments{
\item{alpha}{a vector with two components giving a univariate upper tail probabiliity, and a bivariate upper tail probability, for wwhich the corresponding quanties are desired.}

\item{klist}{a list with compoents k1, k2, k3, and k4, representing first through fourth comulant matrices respectively, having 1 through four dimensions, each dimension 2, respectively.}

\item{nn}{the sample size.  If nn is not one, it is presumed that the cumulants of the sampling distribution under investigation behave in teh usual way; that is, third cumulants are k3/sqrt(nn), and fourth cumulants are k4/(nn.}

\item{xv}{As an alternative to alpha, one may enter ordinates associated with bivariate normal tail probabilities desired, and the tail probabilities are recalcuated from these.}

\item{prestd}{Logical variable indicating whether third and fourth order cumulants are already standardized to unit marginal variance.}
}
\value{
Approximate quantiles.
}
\description{
Calculate the Cornish-Fisher two-dimensional quantile approximation
}
