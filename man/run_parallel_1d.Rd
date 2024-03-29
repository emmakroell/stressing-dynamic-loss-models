% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary_stats.R
\name{run_parallel_1d}
\alias{run_parallel_1d}
\title{Helper for summary_stats()
Runs model at different starting points t0, x0, in parallel}
\usage{
run_parallel_1d(
  niter,
  delta_t,
  t0_vec,
  x0_vec,
  kappa,
  jump_dist,
  stress_type = stress_type,
  stress_parms = stress_parms,
  dt = dt,
  Npaths = Npaths
)
}
\arguments{
\item{niter}{int, number of iterations to take}

\item{delta_t}{float, time increment at which to calculate stats}

\item{t0_vec}{vector, times at which to start model}

\item{x0_vec}{vector, state values at which to start model}

\item{kappa}{double, compound Poisson parameter}

\item{jump_dist}{RPS_dist object, contains jump size distribution}

\item{stress_type}{one of "VaR", "VaR & CVaR"}

\item{stress_parms}{list of parameters for stress
VaR: q: stressed VaR value, c: VaR level, VaR_stress: multiplier for P-VaR
CVAR: q: stressed VaR value, c: VaR level, s: stressed CVaR level,
VaR_stress: multiplier for P-VaR, CVaR_stress: multiplier for P-CVaR
time_stress: time at which stress occurs; if NULL, assumed to be end_time
time steps of conditioning_path}

\item{dt}{float, step size in time}

\item{Npaths}{integer, number of paths}
}
\value{

}
\description{
Helper for summary_stats()
Runs model at different starting points t0, x0, in parallel
}
