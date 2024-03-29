% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/1d_plot_paths.R
\name{plot_paths}
\alias{plot_paths}
\title{plot a RPS object
creates two plots: the first shows elect paths under the stressed model,
along with the chosen quantiles under P and Q, the second shows the
corresponding intensity paths udner Q}
\usage{
plot_paths(
  object,
  Npaths,
  quantiles = list(lower = 0.1, upper = 0.9),
  indices = "random",
  base_size = 14
)
}
\arguments{
\item{object}{RPS_model object}

\item{Npaths}{int, number of paths to plot}

\item{quantiles}{list containing lower and upper, specifies which quantiles
to highlight}

\item{indices}{vector of indices to plot of length Npaths,
or ="random" to choose randomly}
}
\value{

}
\description{
plot a RPS object
creates two plots: the first shows elect paths under the stressed model,
along with the chosen quantiles under P and Q, the second shows the
corresponding intensity paths udner Q
}
