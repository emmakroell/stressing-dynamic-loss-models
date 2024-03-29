% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_parallel.R
\name{run_parallel}
\alias{run_parallel}
\title{Run bivariate VaR stress simulations in parallel for a given set of stresses}
\usage{
run_parallel(stress_seq, alpha, beta, jump_dist, kappa = 5, dt = 0.01)
}
\arguments{
\item{stress_seq}{sequence of stresses}

\item{alpha}{dbl, between 0 and 1 (strictly), alpha level for VaR at time 0.5}

\item{beta}{dbl, between 0 and 1 (strictly), alpha level for VaR at time 1}

\item{jump_dist}{RPS_dist_biv object}

\item{kappa}{dbl, baseline intensity}

\item{dt}{dbl, step size}
}
\value{
tibble of stresses to VaR^P_alpha(X^1_{T/2}) and
resulting VaR^Q_beta(X^1_T + X^2_T), where T=1
}
\description{
Run bivariate VaR stress simulations in parallel for a given set of stresses
}
