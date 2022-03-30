
#' Plot dist of G marginals at multiple points in time and X
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
plot_G_marginals <- function(x_vec,t_vec,kappa,jump_dist,stress_type,stress_parms,N_out=1e4){

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

  # N draws from jump size distribution
  withr::with_seed(720,Draws <- copula::rMvdc(1e5,jump_dist$biv_dist))

  dat <- NULL
  for (t in t_vec){
    res <- sim_G_kappa(t=t,x_vec,eta=eta,kappa=kappa,stress_type=stress_type,
                       stress_parms=stress_parms,dist=jump_dist,Draws=Draws,
                       N_out=N_out)

  # arrange data
    Y1 <- dplyr::as_tibble(t(res$G_Q_1))
    colnames(Y1) <- x_vec
    Y1 <- Y1 %>%
      tidyr::pivot_longer(cols=dplyr::everything(),
                          names_to="x",values_to="value") %>%
      dplyr::mutate(dim=1,time = as.factor(paste("time =",t)))

    Y2 <- dplyr::as_tibble(t(res$G_Q_2))
    colnames(Y2) <- x_vec
    Y2 <- Y2 %>%
      tidyr::pivot_longer(cols=dplyr::everything(),names_to="x",values_to="value") %>%
      dplyr::mutate(dim=2,time = as.factor(paste("time =",t)))

    dat <- rbind(dat,Y1,Y2)
  }

  # plots
  plot_y1 <- dat %>%
    dplyr::filter(dim == 1) %>%
    dplyr::mutate(X=paste("x =",x)) %>%
    ggplot2::ggplot() +
    ggplot2::geom_histogram(ggplot2::aes(x=value,y=..density..,colour=X,fill=X),
                            bins=25,alpha=0.6) +
    ggplot2::geom_line(data=dplyr::as_tibble(cbind(plot_x,plot_y)),
                       ggplot2::aes(x=plot_x,y=plot_y)) +
    ggplot2::facet_grid(X ~ time) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::xlab('') + ggplot2::ylab('') + ggplot2::xlim(c(0,12))

  plot_y2 <- dat %>%
    dplyr::filter(dim == 2) %>%
    dplyr::mutate(X=paste("x =",x)) %>%
    ggplot2::ggplot() +
    ggplot2::geom_histogram(ggplot2::aes(x=value,y=..density..,colour=X,fill=X),
                            bins=25,alpha=0.6) +
    ggplot2::geom_line(data=dplyr::as_tibble(cbind(plot_x,plot_y)),
                       ggplot2::aes(x=plot_x,y=plot_y)) +
    ggplot2::facet_grid(X ~ time) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::xlab('') + ggplot2::ylab('') + ggplot2::xlim(c(0,12))

  return(list(plot_data = dat,
              plot_y1=plot_y1,
              plot_y2=plot_y2))
}


# plot X_t at given time
#' creates two-panel plot of histograms
#'
#' @param object RPS_model object
#' @param time dbl, must be in object$time_vec
#'
#' @return
#' @export
#'
plot_X_histogram <- function(object,time){

  time_index <- match(time,object$time_vec)

  X1 <- object$paths$X1[time_index,]
  X2 <- object$paths$X2[time_index,]

  dplyr::as_tibble(cbind(X1,X2)) %>%
    tidyr::pivot_longer(cols=dplyr::everything(),names_to="X",values_to="value") %>%
    ggplot2::ggplot() +
    ggplot2::geom_histogram(ggplot2::aes(x=value,y=..density..,colour=X,fill=X),
                            bins=20,alpha=0.6) +
    ggplot2::facet_grid(~X) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(legend.position = "none") + ggplot2::xlim(c(0,40))
}


plot_X_point_marginal <- function(object,time){

  time_index <- match(time,object$time_vec)

  plot <- data.frame(X1 = object$paths$X1[time_index,],
                     X2 = object$paths$X2[time_index,]) %>%
    ggplot2::ggplot(aes(X1, X2)) +
    ggplot2::geom_point() +
    ggplot2::theme_bw()

  ggExtra::ggMarginal(plot, type = "histogram")
  # add colour: xparams=list(fill = "orange"), yparams = list(fill="green")
}


# compare X to baseline
#' creates two-panel plot of histograms
#'
#' @param object RPS_model object
#' @param time dbl, must be in object$time_vec
#'
#' @return
#' @export
#'
compare_X_baseline_histogram <- function(object,time){

  time_index <- match(time,object$time_vec)

  X1 <- object$paths$X1[time_index,]
  X2 <- object$paths$X2[time_index,]
  X_baseline <- sim_baseline(object)[time_index,]

  dplyr::as_tibble(cbind(X1,X_baseline,X2)) %>%
    tidyr::pivot_longer(cols=dplyr::everything(),names_to="X",values_to="value") %>%
    ggplot2::ggplot() +
    ggplot2::geom_histogram(ggplot2::aes(x=value,y=..density..,colour=X,fill=X),
                            bins=20,alpha=0.6) +
    ggplot2::geom_vline(ggplot2::aes(xintercept=object$stress_parms$q),lty=2,col="grey")+
    ggplot2::facet_grid(~X) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::xlim(c(0,40))
}


