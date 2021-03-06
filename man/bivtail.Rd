% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bivtail.R
\name{bivtail}
\alias{bivtail}
\title{Calculate the Edgworth two-dimensional tail probability approximation}
\usage{
bivtail(xv, klist, nn = 1, alreadystand = F)
}
\arguments{
\item{xv}{a vector with two components giving the point at which to evaluate the tail approximaton}

\item{klist}{a list with compoents k1, k2, k3, and k4, representing first through fourth comulant matrices respectively, having 1 through four dimensions, each dimension 2, respectively.}

\item{nn}{the sample size.  If nn is not one, it is presumed that the cumulants of the sampling distribution under investigation behave in teh usual way; that is, third cumulants are k3/sqrt(nn), and fourth cumulants are k4/(nn.}

\item{alreadystand}{Logical flag, false if ordinates need to be standardized to zero mean and unit variance, and true otherwise.}
}
\value{
A two-component vector with the marginal univariate Edgeworth approximation to the tail for the first dimension, and the bivariate tail approximation.
}
\description{
Calculate the Edgworth two-dimensional tail probability approximation
}
