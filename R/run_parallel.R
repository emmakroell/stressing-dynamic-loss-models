#' Run bivariate VaR stress simulations in parallel for a given set of stresses
#'
#' @param stress_seq sequence of stresses
#' @param alpha dbl, between 0 and 1 (strictly), alpha level for VaR at time 0.5
#' @param beta dbl, between 0 and 1 (strictly), alpha level for VaR at time 1
#' @param jump_dist RPS_dist_biv object
#' @param kappa dbl, baseline intensity
#' @param dt dbl, step size
#'
#' @return tibble of stresses to VaR^P_alpha(X^1_{T/2}) and
#'           resulting VaR^Q_beta(X^1_T + X^2_T), where T=1
#' @export
#'
run_parallel <- function(stress_seq,alpha,beta,jump_dist,kappa=5,dt=1e-2) {
  tm <- proc.time()

  # Start parallel computing
  mc.cores <- getOption("mc.cores", max(1, detectCores()-4))
  cl <- makeCluster(mc.cores,outfile = "")
  registerDoParallel(cl)

  res <- foreach(stress=stress_seq,
                 .export = c("stressed_sim_biv","sim_G_kappa_biv",
                             "compute_VaR", "eta_VaR",
                             "FST","h_VaR","P_prob","RPS_model"),
                 .combine='rbind') %dopar% {
                   sim <- stressed_sim_biv(kappa = kappa,
                                           jump_dist = jump_dist,
                                           stress_type = "VaR",
                                           stress_parms = list(c=alpha,
                                                               VaR_stress=stress,
                                                               time_stress=0.5),
                                           Npaths=1e4, end_time=1, dt=dt, beep=FALSE)

                   X1_VaR <- quantile(sim$paths$X1[length(sim$time_vec),],
                                      beta,type=1)
                   X2_VaR <- quantile(sim$paths$X2[length(sim$time_vec),],
                                      beta,type=1)
                   sum_VaR <- quantile(sim$paths$X1[length(sim$time_vec),] +
                                         sim$paths$X2[length(sim$time_vec),],
                                       beta,type=1)

                   return(dplyr::as_tibble(list(alpha = alpha,
                                 stress = stress,
                                 X1_VaR = as.numeric(X1_VaR),
                                 X2_VaR = as.numeric(X2_VaR),
                                 sum_VaR = as.numeric(sum_VaR))))
                 }

  stopCluster(cl)
  cat("Run time:", (proc.time() - tm)[3]/60, "min \n")
  beepr::beep()
  return(res)
}


#' Find optimal stress to VaR^P_alpha(X^1_{T/2}) to achieve a set stress to
#'  VaR^Q_beta(X^1_T + X^2_T)
#'
#' @param floor dbl, floor of interval of stresses to search over, must be
#'                   below true stress
#' @param ceiling dbl, ceiling of interval of stresses to search over, must be
#'                    above true stress
#' @param perc_target target stress to VaR^Q_beta(X^1_T + X^2_T), as a percentage
#'                    i.e. 10, not 0.1
#' @param sum_VaR_P   VaR_beta(X^1_T + X^2_T) under the baseline measure
#' @param alpha dbl, between 0 and 1 (strictly), alpha level for VaR at time 0.5
#' @param beta dbl, between 0 and 1 (strictly), alpha level for VaR at time 1
#' @param jump_dist RPS_dist_biv object
#' @param kappa dbl, baseline intensity
#' @param dt dbl, step size
#' @param tol dbl, tolerance at which to stop search
#' @param max_iter int, max number of iteration at which to stop search
#'
#' @return $table: tibble of stresses to VaR^P_alpha(X^1_{T/2}) and
#'           resulting VaR^Q_beta(X^1_T + X^2_T), where T=1
#'         $stress: optimal stress
#' @export
#'
find_stress_parallel <- function(floor,ceiling,perc_target,sum_VaR_P,
                             alpha,beta,jump_dist,kappa=5,dt=1e-2,
                             tol=1e-6,max_iter=4) {
  # initialize
  output <- NULL
  iter <- 0
  diff <- 1
  ncores <- max(1, detectCores()-4) # 6 cores

  while (ceiling - floor > tol & iter < max_iter){
    stress_seq <- seq(floor,ceiling,length.out=ncores)

    # compute VaR on grid
    res <- run_parallel(stress_seq=stress_seq, alpha=alpha,
                        beta=beta, jump_dist=jump_dist,
                        kappa=kappa, dt=dt)

    res_mod <- res %>% dplyr::mutate(sum_VaR_P = sum_VaR_P) %>%
      dplyr::mutate(perc_incr = (sum_VaR - sum_VaR_P)/sum_VaR_P * 100) %>%
      arrange(perc_incr)

    # find floor and ceiling
    floor_ind <- tail(which(res_mod$perc_incr < perc_target),1)
    ceiling_ind <- head(which(res_mod$perc_incr >= perc_target),1)

    # return error if interval is wrong
    if (length(floor_ind) == 0) {
      stop("Solution not in given interval.")
    }
    if (floor_ind == length(stress_seq)) {
      stop("Solution not in given interval.")
    }

    # update
    floor <- res_mod$stress[floor_ind]
    ceiling <- res_mod$stress[ceiling_ind]
    output <- rbind(output,res_mod)
    iter <- iter + 1

    cat("iter = ", iter, ", floor = ", floor, ", ceiling ", ceiling,
        "\n f(floor) = ",res_mod$perc_incr[floor_ind],
        "f(ceiling) = ",res_mod$perc_incr[ceiling_ind], "\n")

    if (floor > ceiling) stop("Floor greater than ceiling.")
  }

  # return midpoint
  stress <- (floor + ceiling) / 2

  return(list(table = output,
              stress = stress,
              niter = iter,
              x_diff = ceiling - floor,
              fx_diff = (res_mod$perc_incr[floor_ind] +
                           res_mod$perc_incr[ceiling_ind])/2 - perc_target))
}
