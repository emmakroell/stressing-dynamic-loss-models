% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class.R
\name{RPS_dist_univ}
\alias{RPS_dist_univ}
\title{Defines RPS_dist class, class giving a distribution}
\usage{
RPS_dist_univ(dist_fun, dens_fun, sim_fun, char_fun, mean_fun, min, max, parms)
}
\arguments{
\item{dist_fun}{function of x and parms, returns distribution function}

\item{dens_fun}{function of x and parms, returns density function}

\item{sim_fun}{function of x and parms, returns random draw generation function}

\item{char_fun}{function of x and parms, returns characteristic function}

\item{mean_fun}{function of parms, returns mean}

\item{min}{dbl, where to start integration of random variable}

\item{max}{dbl, where to end integration of random variable}

\item{parms}{list, parameters for distribution}
}
\value{

}
\description{
Defines RPS_dist class, class giving a distribution
}
