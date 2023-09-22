

#' Helper for summary_stats()
#' Runs model at different starting points t0, x0, in parallel
#'
#' @param niter int, number of iterations to take
#' @param delta_t float, time increment at which to calculate stats
#' @param t0_vec vector, times at which to start model
#' @param x0_vec vector, state values at which to start model
#' @param kappa double, compound Poisson parameter
#' @param jump_dist RPS_dist object, contains jump size distribution
#' @param stress_type one of "VaR", "VaR & CVaR"
#' @param stress_parms list of parameters for stress
#' VaR: q: stressed VaR value, c: VaR level, VaR_stress: multiplier for P-VaR
#' CVAR: q: stressed VaR value, c: VaR level, s: stressed CVaR level,
#' VaR_stress: multiplier for P-VaR, CVaR_stress: multiplier for P-CVaR
#' time_stress: time at which stress occurs; if NULL, assumed to be end_time
#'        time steps of conditioning_path
#' @param dt float, step size in time
#' @param Npaths integer, number of paths
#'
#' @return
#' @export
#'
#' @examples
run_parallel_1d <- function(niter, delta_t,
                            t0_vec, x0_vec,
                            kappa, jump_dist,
                            stress_type = stress_type,
                            stress_parms = stress_parms,
                            dt=dt, Npaths = Npaths) {
  tm <- proc.time()


  # Start parallel computing
  mc.cores <- getOption("mc.cores", max(1, detectCores()-4))
  cl <- makeCluster(mc.cores,outfile = "")
  registerDoParallel(cl)

  res <- foreach(i=1:niter,
                 .export = c("stressed_sim", "sim_baseline","sim_G_kappa",
                             "compute_CVaR", "compute_VaR", "CVaR_calc",
                             "eta_CVaR", "eta2_eq", "P_expectation",
                             "FST","h_CVaR","P_prob","RPS_model"),
                 .combine='rbind') %dopar% {
                   sim <- stressed_sim(kappa = kappa, jump_dist = jump_dist,
                                       stress_type = stress_type,
                                       stress_parms = stress_parms,
                                       X0 = x0_vec[i], t0 = t0_vec[i],
                                       end_time=t0_vec[i] + delta_t,
                                       dt=dt, Npaths = Npaths, beep = FALSE)

                   results_stressed <- sim$paths[dim(sim$paths)[1],] - x0_vec[i]

                   # baseline
                   baseline_sim <- sim_baseline(sim)
                   results_baseline <- baseline_sim[dim(baseline_sim)[1],] - x0_vec[i]

                  # compute the risk measures of both the stressed and baseline runs
                   return(dplyr::as_tibble(list(stressed_median = quantile(results_stressed, probs=0.5,type=1),
                                                baseline_median = quantile(results_baseline, probs=0.5,type=1),
                                                stressed_mean = mean(results_stressed),
                                                baseline_mean = mean(results_baseline),
                                                stressed_VaR = quantile(results_stressed, probs=0.9,type=1),
                                                baseline_VaR = quantile(results_baseline, probs=0.9,type=1),
                                                stressed_CVaR = CVaR_calc(results_stressed, alpha=0.9),
                                                baseline_CVaR = CVaR_calc(results_baseline, alpha=0.9),
                                                x0 = rep(x0_vec[i]),
                                                t0 = rep(t0_vec[i]))))
                 }

  stopCluster(cl)
  cat("Run time:", (proc.time() - tm)[3]/60, "min \n")
  beepr::beep()
  return(res)
}

#' Compute risk measures of X_{t+delta} - X_{t-} given X_{t-}
#' for a grid of times t
#'
#' @param t_grid float, starting time
#' @param delta_t float, time increment at which to calculate stats
#' @param kappa double, compound Poisson parameter
#' @param jump_dist RPS_dist object, contains jump size distribution
#' @param stress_type one of "VaR", "VaR & CVaR"
#' @param stress_parms list of parameters for stress
#' VaR: q: stressed VaR value, c: VaR level, VaR_stress: multiplier for P-VaR
#' CVAR: q: stressed VaR value, c: VaR level, s: stressed CVaR level,
#' VaR_stress: multiplier for P-VaR, CVaR_stress: multiplier for P-CVaR
#' time_stress: time at which stress occurs; if NULL, assumed to be end_time
#' @param conditioning_path vector, path to condition on (run of the model)
#' @param conditioning_path_times vector, time vector that matches the
#'        time steps of conditioning_path
#' @param dt float, step size in time
#' @param Npaths integer, number of paths
#'
#'
#' @return
#' @export
#'
summary_stats <- function(t_grid, delta_t, kappa, jump_dist, stress_type,
                          stress_parms, conditioning_path,
                          conditioning_path_times, dt=1e-2, Npaths=1e4){

  # times at which to evaluate
  t_eval_vec <- seq(0,1-delta_t,by=t_grid)

  # check that eval points
  stopifnot(t_grid/dt == round(t_grid/dt))
  # find correct points in the conditioning paths that match t_eval_vec
  eval_points <- which(round(conditioning_path_times,4) %in% round(t_eval_vec,4))

  # simulate model at t0, x0 in parallel
  res <- run_parallel_1d(niter = length(eval_points),
                         delta_t = delta_t,
                         t0_vec = t_eval_vec,
                         x0_vec = conditioning_path[eval_points],
                         kappa = kappa, jump_dist = jump_dist,
                         stress_type = stress_type,
                         stress_parms = stress_parms,
                         dt=dt, Npaths = Npaths)
  return(res)
}


#
#' helper function for CVaR
#'
#' @param X vector, data
#' @param alpha float, CVaR level
#'
#' @return
#' @export
#'
CVaR_calc <- function(X,alpha){
  VaR <- quantile(X,alpha,type=1)
  mean(X[X>VaR])
}
