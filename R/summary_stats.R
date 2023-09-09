
#' Compute summary statistics at X_{t_delta} given X_{t0} = X0
#'
#' @param X0 vector, initial conditions for X
#' @param t0 float, starting time
#' @param t_delta float, time at which to calculate stats
#' @param kappa double, compound Poisson parameter
#' @param jump_dist RPS_dist object, contains jump size distribution
#' @param stress_type one of "VaR", "VaR & CVaR"
#' @param stress_parms list of parameters for stress
#' VaR: q: stressed VaR value, c: VaR level, VaR_stress: multiplier for P-VaR
#' CVAR: q: stressed VaR value, c: VaR level, s: stressed CVaR level,
#' VaR_stress: multiplier for P-VaR, CVaR_stress: multiplier for P-CVaR
#' time_stress: time at which stress occurs; if NULL, assumed to be end_time
#' @param end_time float, when to end sim
#' @param dt float, step size in time
#' @param Npaths integer, number of paths
#'
#' @return
#' @export
#'
summary_stats <- function(X0, t0, t_delta, kappa, jump_dist, stress_type,
                     stress_parms, end_time=1, dt=1e-2, Npaths=1e4){

  # delete results_all
  results_all <- NULL

  # save results
  results_condensed <- NULL

  for (i in 1:length(X0)){
    x <- X0[i]

    # run model
    stressed_res <- stressed_sim(kappa = kappa, jump_dist = jump_dist,
                             stress_type = stress_type,
                             stress_parms = stress_parms,
                             X0 = x, t0 = t0,
                             end_time=end_time, dt=dt, Npaths = Npaths)

    # pull out draws from stressed & basline models at time t_delta
    stressed <- stressed_res$paths[which(stressed_res$time_vec == t_delta),]
    baseline <- sim_baseline(stressed_res)[which(stressed_res$time_vec == t_delta),]

    results_all <- cbind(results_all,stressed,baseline)

    # calculate summary statistics
    res <- cbind(stressed,baseline) |>
      dplyr::as_tibble() |>
      dplyr::summarise(dplyr::across(dplyr::everything(),
                                     list(mean = mean,
                                          VaR  = ~ quantile(.x, probs=0.9,type=1),
                                          CVaR  = ~ CVaR_calc(.x,alpha=0.9)
                                     ))) |>
      cbind(x) |>
      dplyr::as_tibble()

    # update results
    results_condensed <- rbind(results_condensed,res)

  }

  # pivot longer
  results_condensed <- results_condensed |>
    tidyr::pivot_longer(cols = 1:6,
                        names_to = c(".value", "statistic"),
                        names_sep="_")

  return(list(results_all=results_all,
              results_condensed=results_condensed))
}

# helped func for CVaR
CVaR_calc <- function(X,alpha){
  VaR <- quantile(X,alpha,type=1)
  mean(X[X>VaR])
}

