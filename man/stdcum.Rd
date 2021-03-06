% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stdcum.R
\name{stdcum}
\alias{stdcum}
\title{Standardize a list of cumulants to have marginal zero mean and unit variance.}
\usage{
stdcum(klist)
}
\arguments{
\item{klist}{a list with compoents k1, k2, k3, and k4, representing first through fourth comulant matrices respectively, having 1 through four dimensions, each dimension 2, respectively.}
}
\value{
A list with the same components, with zero mean, unit marginal variance, and adjusted cumulants.  No adjustment is made for correlation.
}
\description{
Standardize a list of cumulants to have marginal zero mean and unit variance.
}
