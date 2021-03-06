% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/givex.R
\name{fun.givex}
\alias{fun.givex}
\title{Give the bivariate normal upper quantile}
\usage{
fun.givex(alpha, rho, expect = NULL, sd = NULL)
}
\arguments{
\item{alpha}{a vector with two components giving upper quantiles}

\item{rho}{the correlation for the bivariate normal with two normal components.}

\item{expect}{the two-component vector of expected values of the bivariate distribution, if not zero.}

\item{sd}{the two-component vector of marginal standard deviations of the bivariate distribution, if not 1.}
}
\value{
A two-component vector whose first component is the upper standard normal quantile associated with alpha[1], and whose second component is that one to make the bivariate normal upper tail area alpha[2].
}
\description{
Give the bivariate normal upper quantile
}
