#' Plot kappa^Q in 3D using plotly
#'
#' @param kappa double, compound Poisson parameter
#' @param jump_dist RPS_dist object, contains jump size distribution
#' @param stress_type one of "VaR", "VaR & CVaR"
#' @param stress_parms list of parameters for stress
#' VaR: q: stressed VaR value, c: VaR level, VaR_stress: multiplier for P-VaR
#' CVAR: q: stressed VaR value, c: VaR level, s: stressed CVaR level,
#' VaR_stress: multiplier for P-VaR, CVaR_stress: multiplier for P-CVaR
#' @param endtime float, when to end sim
#' @param dt float, step size in time
#'
#' @return
#' @export
#'
plot_kappa <- function(kappa,jump_dist,stress_type,stress_parms,endtime=1,dt=1e-2){

  if (is.null(stress_parms$q) & !(is.null(stress_parms$VaR_stress))){
    stress_parms$q <- stress_parms$VaR_stress * compute_VaR(kappa,stress_parms$c,jump_dist)
  }
  if ((stress_type == "CVaR") & is.null(stress_parms$s) &
      !(is.null(stress_parms$CVaR_stress))){
    stress_parms$s <- stress_parms$CVaR_stress *
      compute_CVaR(kappa=kappa,q=stress_parms$q,c=stress_parms$c,dist=jump_dist)
    if (stress_parms$s < stress_parms$q) {
      stop("Incompatible parameter choice: attempt to set CVaR below VaR.\n")
    }
  }

  eta <- switch(stress_type,
                "VaR" = eta_VaR(kappa=kappa, stress_parms=stress_parms,
                                dist=jump_dist),
                "CVaR" = eta_CVaR(kappa=kappa, stress_parms=stress_parms,
                                  dist=jump_dist))

  set.seed(720)
  Draws <- jump_dist$sim_fun(1e4,jump_dist$parms)  # N draws

  Nsteps <- endtime / dt + 1   # number of steps to take
  # create grid of kappa
  times <- seq(0,endtime,by=dt)
  # choose X max based on q:
  x.max <- stress_parms$q + 5
  x.seq <- seq(0,x.max,by=0.1)
  len.x <- length(x.seq)

  kappa_Q <- matrix(nrow = Nsteps, ncol = len.x)
  for (i in 1:Nsteps){
    kappa_Q[i,] <- sim_G_kappa(x=x.seq,t=times[i],eta=eta,kappa=kappa,
                               stress_type = stress_type,
                               stress_parms=stress_parms,
                               dist=jump_dist,Draws=Draws)$kappa_Q
  }

  # make the plot
  axx <- list(
    ticketmode = 'array',
    ticktext = seq(0,x.max,by=5),
    tickvals = 10*seq(0,x.max,by=5),
    range = 10*c(0,x.max)
  )
  # y-axis ticks
  axy <- list(
    ticketmode = 'array',
    ticktext = seq(0,endtime,by=0.25),
    tickvals = (1/dt)*seq(0,endtime,by=0.25),
    range = (1/dt)*c(0,endtime),
    title = "time"
  )
  # z-axis ticks
  axz <- list(
    nticks = 6,
    range = c(min(kappa_Q),max(kappa_Q)),
    title = ""
  )

  # make plot
  fig <- plotly::plot_ly(z = kappa_Q, colorscale ='Jet') %>%
    plotly::add_surface() %>%
    plotly::layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz))

  return(list(kappa_Q=kappa_Q,fig=fig))
}


# res <- plot_kappa(kappa=5,jump_dist=gamma_2_1,stress_type="CVaR",
#            stress_parms = list(c=0.9,VaR_stress=1.1,CVaR_stress=1.08))
# res <- plot_kappa(kappa=5,jump_dist=gamma_2_1,stress_type="VaR",
#                   stress_parms = list(c=0.9,VaR_stress=1.1))

#' Plot dist of G^Q at multiple points in time and X
#'
#' @param x_vec vector of state variable values
#' @param t_vec vector of dbl, between 0 and 1
#' @param kappa double, compound Poisson parameter
#' @param jump_dist RPS_dist object, contains jump size distribution
#' @param stress_type one of "VaR", "VaR & CVaR"
#' @param stress_parms list of parameters for stress
#' VaR: q: stressed VaR value, c: VaR level, VaR_stress: multiplier for P-VaR
#' CVAR: q: stressed VaR value, c: VaR level, s: stressed CVaR level,
#' VaR_stress: multiplier for P-VaR, CVaR_stress: multiplier for P-CVaR
#' @param N_out int
#'
#' @return
#' @export
#'
plot_G_Q <- function(x_vec,t_vec,kappa,jump_dist,stress_type,stress_parms,N_out=1e4){

  plot_x <- seq(0,12,0.1)
  plot_y <- jump_dist$dens_fun(plot_x,jump_dist$parms)

  if (is.null(stress_parms$q) & !(is.null(stress_parms$VaR_stress))){
    stress_parms$q <- stress_parms$VaR_stress * compute_VaR(kappa,stress_parms$c,jump_dist)
  }
  if ((stress_type == "CVaR") & is.null(stress_parms$s) &
      !(is.null(stress_parms$CVaR_stress))){
    stress_parms$s <- stress_parms$CVaR_stress *
      compute_CVaR(kappa=kappa,q=stress_parms$q,c=stress_parms$c,dist=jump_dist)
    if (stress_parms$s < stress_parms$q) {
      stop("Incompatible parameter choice: attempt to set CVaR below VaR.\n")
    }
  }

  eta <- switch(stress_type,
                "VaR" = eta_VaR(kappa=kappa, stress_parms=stress_parms,
                                dist=jump_dist),
                "CVaR" = eta_CVaR(kappa=kappa, stress_parms=stress_parms,
                                  dist=jump_dist))

  set.seed(720)
  Draws <- jump_dist$sim_fun(1e4,jump_dist$parms)  # N draws


  dat <- NULL

  for (t in t_vec){
    res <- sapply(x_vec,
                  function(x) sim_G_kappa(t=t,x,eta=eta,kappa=kappa,
                                          stress_type=stress_type,
                                          stress_parms=stress_parms,
                                          dist=jump_dist,Draws=Draws,
                                          N_out=N_out)$G_Q)

    dat <- rbind(dat,dplyr::mutate(dplyr::as_tibble(res),
                                   time = as.factor(paste("time =",t))))
  }

  colnames(dat) <- c(x_vec,"time")

  dat %>%
    tidyr::pivot_longer(cols = !c("time"),names_to = "X",values_to = "value") %>%
    dplyr::mutate(X = factor(X, levels = x_vec)) %>%
    dplyr::mutate(X=paste("X =",X)) %>%
    ggplot2::ggplot() +
    ggplot2::geom_histogram(ggplot2::aes(x=value,y=..density..,colour=X,fill=X),
                            bins=25,alpha=0.6) +
    ggplot2::geom_line(data=dplyr::as_tibble(cbind(plot_x,plot_y)),
                       ggplot2::aes(x=plot_x,y=plot_y)) +
    ggplot2::facet_grid(X ~ time) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(legend.position = "none") + ggplot2::xlab('') + ggplot2::ylab('')

}


# plot_G_Q(x_vec=c(15,17,19,21),t=c(0.5,1),kappa=5,jump_dist=gamma_2_1,
#          stress_type="VaR",stress_parms = list(c=0.9,VaR_stress=1.1))
# plot_G_Q(x_vec=c(15,17,19,21),t=c(0.5,1),kappa=5,jump_dist=gamma_2_1,
#          stress_type="CVaR",stress_parms = list(c=0.9,VaR_stress=1.1,CVaR_stress=1.08))


