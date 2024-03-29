% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FST_implementations.R
\name{eta_VaR}
\alias{eta_VaR}
\title{Computes the optimal eta for VaR stress using FST}
\usage{
eta_VaR(kappa, stress_parms, dist)
}
\arguments{
\item{kappa}{double, compound Poisson parameter}

\item{stress_parms}{list of parameters for stress; must contain:
q: stressed VaR value, time_stress: time at which stress occurs}

\item{dist}{RPS_dist object, contains jump size distribution}
}
\value{
double, optimal eta for VaR stress
}
\description{
Computes the optimal eta for VaR stress using FST
}
