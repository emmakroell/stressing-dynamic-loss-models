# Mixture model code
new_RPS_dist_mix <- function(sim_fun_a, sim_fun_b, dens_fun_a, dens_fun_b,
                                 char_fun, mean_fun,parms) {

  model <- list(sim_fun_a = sim_fun_a,
                sim_fun_b = sim_fun_b,
                dens_fun_a = dens_fun_a,
                dens_fun_b = dens_fun_b,
                char_fun = char_fun,
                mean_fun = mean_fun,
                parms = parms
  )
  ## Name of the class
  attr(model, "class") <- "RPS_dist_mixture"
  return(model)
}
# need to include p (weight) in parms

mean.RPS_dist_mixture <- function(object) {
  with(object, mean_fun(parms))
}

new_RPS_model_mix <- function(jump_dist, kappa, stress_type, stress_parms,
                          paths, time_vec, kappa_Q, jumps, p_Q){

  stress_type <- match.arg(stress_type, c("VaR", "CVaR"))

  model <- list(jump_dist = jump_dist,
                kappa = kappa,
                stress_type = stress_type, # type of stress: VaR or VaR & CVaR
                stress_parms = stress_parms, # list of stress parameters corresponding to stress type
                paths = paths,
                time_vec = time_vec,
                kappa_Q = kappa_Q,
                p_Q = p_Q,
                jumps = jumps
  )
  ## Name of the class
  attr(model, "class") <- "RPS_model"
  return(model)
}


#' Simulate the stressed model
#'
#' @param kappa double, compound Poisson parameter
#' @param jump_dist RPS_dist object, contains jump size distribution
#' @param stress_type one of "VaR", "VaR & CVaR"
#' @param stress_parms list of parameters for stress
#' VaR: q: stressed VaR value, c: VaR level, VaR_stress: multiplier for P-VaR
#' CVAR: q: stressed VaR value, c: VaR level, s: stressed CVaR level,
#' VaR_stress: multiplier for P-VaR, CVaR_stress: multiplier for P-CVaR
#' @param Npaths integer, number of paths
#' @param endtime float, when to end sim
#' @param dt float, step size in time
#'
#' @return RPS_model object
#' @export
stressed_sim_mix <- function(kappa, jump_dist, stress_type = "VaR",
                             stress_parms = list(c=c,q=q,s=s,
                                                 VaR_stress=VaR_stress,
                                                 CVaR_stress=CVaR_stress),
                             Npaths=1e4,endtime=1, dt=1e-2){

  # first, draw from the jump size distribution
  Ndraws <- 1e4
  with(jump_dist,{
    withr::with_seed(720,{
      # N draws from both jump size distributions
      Draws_a <- sim_fun_a(Ndraws,parms)
      Draws_b <- sim_fun_b(Ndraws,parms)

      Nsteps <- endtime / dt + 1   # number of steps to take
      times <- seq(0,endtime,by=dt)

      if (is.null(stress_parms$q) & !(is.null(stress_parms$VaR_stress))){
        # baseline <- sim_baseline_mixture(jump_dist=jump_dist,kappa=kappa,
        #                                  Npaths=Npaths,time_vec=times)
        # VaR_estimate <- quantile(baseline$X1[length(times),], 0.9, type=1)
        # stress_parms$q <- stress_parms$VaR_stress * VaR_estimate
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

      # # create grid of Gs and kappas
      # # choose X max based on parameters:
      # x_max <- max(6 * kappa * mean(jump_dist),60) # IMPROVE LATER
      # x1_seq <- seq(0,x_max,by=1e-2)
      # x2_seq <- seq(0,x_max,by=1e-2)
      # len_x <- length(x1_seq)
      #
      # G_grid1 <- array(dim=c(Nsteps,len_x))
      # G_grid2 <- array(dim=c(Nsteps,len_x))
      # kappa_grid <- matrix(nrow = Nsteps, ncol = len_x)
      # p_grid <- matrix(nrow = Nsteps, ncol = len_x)
      #
      # # compute G and kappa at grid points
      # for (i in 1:Nsteps){
      #   distorted <- sim_G_kappa_mix(t=times[i],x=x1_seq,eta=eta,kappa=kappa,
      #                                stress_type=stress_type,
      #                                stress_parms=stress_parms,p=parms$p,
      #                                dist=jump_dist,Draws_a=Draws_a,
      #                                Draws_b=Draws_b)
      #   G_grid1[i,] <- distorted$G_Q_1
      #   G_grid2[i,] <- distorted$G_Q_2
      #   kappa_grid[i,] <- distorted$kappa_Q
      #   p_grid[i,] <- distorted$p_Q
      # }

      # initialize
      X1 <- matrix(nrow=Nsteps,ncol=Npaths) # empty matrix to store results
      X2 <- matrix(nrow=Nsteps,ncol=Npaths)
      kappa_Q <- matrix(nrow=Nsteps,ncol=Npaths)
      p_Q <- matrix(nrow=Nsteps,ncol=Npaths)

      jump_counter <- matrix(nrow=Nsteps,ncol=Npaths)
      jump_counter[1,] <- rep(0,Npaths)

      X1[1,] <- rep(0,Npaths)   # X starts at 0 for all paths
      X2[1,] <- rep(0,Npaths)
      times <- seq(0,endtime,by=dt)

      for (i in 2:Nsteps){
        distorted <- sim_G_kappa_mix(t=times[i],x=X1[i-1,],eta=eta,
                                     kappa=kappa,stress_type=stress_type,
                                     stress_parms=stress_parms,p=parms$p,
                                     dist=jump_dist,Draws_a=Draws_a,
                                     Draws_b=Draws_b)
        U <- runif(Npaths)  # generate Npaths independent unif[0,1]

        # compute kappa.Q and draw from G.Q, using interpolation
        kappa_Q[i,] <- distorted$kappa_Q  # store kappa
        p_Q[i,] <- distorted$p_Q # store p
        G_Q_draw1 <- distorted$G_Q_1
        G_Q_draw2 <- distorted$G_Q_2
        # browser()
        # simulate forward:
        # X1[i,] <- X1[i-1,] + G_Q_draw1 * as.integer(U < (1 - exp(-kappa_Q[i,] * dt)))
        # X2[i,] <- X2[i-1,] + G_Q_draw2 * as.integer(U < (1 - exp(-kappa_Q[i,] * dt)))
        jump <- as.integer(U < (1 - exp(-kappa_Q[i,] * dt)))
        X1[i,] <- X1[i-1,] + G_Q_draw1 * jump
        X2[i,] <- X2[i-1,] + G_Q_draw2 * jump
        jump_counter[i,] <- jump_counter[i-1,] + jump

        # # if X passed the grid max, return an error:
        # if (any(c(X1[i,],X2[i,]) > x_max)){
        #    cat("max X1 ", max(X1[i,]), "max X2 ", max(X2[i,]), "\n")
        #    stop("Path value exceeded grid max. Increase x max. \n")
        #   }
      }
    })
    beepr::beep()
    new_RPS_model_mix(jump_dist = jump_dist, kappa = kappa,
                  stress_type = stress_type, stress_parms = stress_parms,
                  time_vec = times, paths = list(X1=X1, X2=X2),
                  kappa_Q = kappa_Q, p_Q = p_Q, jumps = jump_counter)
  })
}


# OLD
# stressed_sim_mix <- function(kappa, jump_dist, stress_type = "VaR",
#                              stress_parms = list(c=c,q=q,s=s,
#                                                  VaR_stress=VaR_stress,
#                                                  CVaR_stress=CVaR_stress),
#                              Npaths=1e4,endtime=1, dt=1e-2){
#
#   # first, draw from the jump size distribution
#   Ndraws <- 1e4
#   with(jump_dist,{
#     withr::with_seed(720,{
#       # N draws from both jump size distributions
#       Draws_a <- sim_fun_a(Ndraws,parms)
#       Draws_b <- sim_fun_b(Ndraws,parms)
#
#       Nsteps <- endtime / dt + 1   # number of steps to take
#       times <- seq(0,endtime,by=dt)
#
#       if (is.null(stress_parms$q) & !(is.null(stress_parms$VaR_stress))){
#         # baseline <- sim_baseline_mixture(jump_dist=jump_dist,kappa=kappa,
#         #                                  Npaths=Npaths,time_vec=times)
#         # VaR_estimate <- quantile(baseline$X1[length(times),], 0.9, type=1)
#         # stress_parms$q <- stress_parms$VaR_stress * VaR_estimate
#         stress_parms$q <- stress_parms$VaR_stress * compute_VaR(kappa,stress_parms$c,jump_dist)
#       }
#       if ((stress_type == "CVaR") & is.null(stress_parms$s) &
#           !(is.null(stress_parms$CVaR_stress))){
#         stress_parms$s <- stress_parms$CVaR_stress *
#           compute_CVaR(kappa=kappa,q=stress_parms$q,c=stress_parms$c,dist=jump_dist)
#         if (stress_parms$s < stress_parms$q) {
#           stop("Incompatible parameter choice: attempt to set CVaR below VaR.\n")
#         }
#       }
#       eta <- switch(stress_type,
#                     "VaR" = eta_VaR(kappa=kappa, stress_parms=stress_parms,
#                                     dist=jump_dist),
#                     "CVaR" = eta_CVaR(kappa=kappa, stress_parms=stress_parms,
#                                       dist=jump_dist))
#
#       # create grid of Gs and kappas
#       # choose X max based on parameters:
#       x_max <- max(6 * kappa * mean(jump_dist),60) # IMPROVE LATER
#       x1_seq <- seq(0,x_max,by=1e-2)
#       x2_seq <- seq(0,x_max,by=1e-2)
#       len_x <- length(x1_seq)
#
#       G_grid1 <- array(dim=c(Nsteps,len_x))
#       G_grid2 <- array(dim=c(Nsteps,len_x))
#       kappa_grid <- matrix(nrow = Nsteps, ncol = len_x)
#       p_grid <- matrix(nrow = Nsteps, ncol = len_x)
#
#       # compute G and kappa at grid points
#       for (i in 1:Nsteps){
#         distorted <- sim_G_kappa_mix(t=times[i],x=x1_seq,eta=eta,kappa=kappa,
#                                      stress_type=stress_type,
#                                      stress_parms=stress_parms,p=parms$p,
#                                      dist=jump_dist,Draws_a=Draws_a,
#                                      Draws_b=Draws_b)
#         G_grid1[i,] <- distorted$G_Q_1
#         G_grid2[i,] <- distorted$G_Q_2
#         kappa_grid[i,] <- distorted$kappa_Q
#         p_grid[i,] <- distorted$p_Q
#       }
#
#       # initialize
#       X1 <- matrix(nrow=Nsteps,ncol=Npaths) # empty matrix to store results
#       X2 <- matrix(nrow=Nsteps,ncol=Npaths)
#       kappa_Q <- matrix(nrow=Nsteps,ncol=Npaths)
#       p_Q <- matrix(nrow=Nsteps,ncol=Npaths)
#
#       jump_counter <- matrix(nrow=Nsteps,ncol=Npaths)
#       jump_counter[1,] <- rep(0,Npaths)
#
#       X1[1,] <- rep(0,Npaths)   # X starts at 0 for all paths
#       X2[1,] <- rep(0,Npaths)
#       times <- seq(0,endtime,by=dt)
#
#       for (i in 2:Nsteps){
#         U <- runif(Npaths)  # generate Npaths independent unif[0,1]
#
#         # compute kappa.Q and draw from G.Q, using interpolation
#         kappa_Q[i,] <- approx(x=x1_seq,y=kappa_grid[i-1,],xout=X1[i-1,])$y  # store kappa
#         p_Q[i,] <- approx(x=x1_seq,y=p_grid[i-1,],xout=X1[i-1,])$y  # store p
#         G_Q_draw1 <- approx(x=x1_seq,y=G_grid1[i-1,],xout=X1[i-1,])$y
#         G_Q_draw2 <- approx(x=x2_seq,y=G_grid2[i-1,],xout=X2[i-1,])$y
#         browser()
#         # simulate forward:
#         # X1[i,] <- X1[i-1,] + G_Q_draw1 * as.integer(U < (1 - exp(-kappa_Q[i,] * dt)))
#         # X2[i,] <- X2[i-1,] + G_Q_draw2 * as.integer(U < (1 - exp(-kappa_Q[i,] * dt)))
#         jump <- as.integer(U < (1 - exp(-kappa_Q[i,] * dt)))
#         X1[i,] <- X1[i-1,] + G_Q_draw1 * jump
#         X2[i,] <- X2[i-1,] + G_Q_draw2 * jump
#         jump_counter[i,] <- jump_counter[i-1,] + jump
#
#         # if X passed the grid max, return an error:
#         if (any(c(X1[i,],X2[i,]) > x_max)){
#           cat("max X1 ", max(X1[i,]), "max X2 ", max(X2[i,]), "\n")
#           stop("Path value exceeded grid max. Increase x max. \n")
#         }
#       }
#     })
#     beepr::beep()
#     new_RPS_model_mix(jump_dist = jump_dist, kappa = kappa,
#                       stress_type = stress_type, stress_parms = stress_parms,
#                       time_vec = times, paths = list(X1=X1, X2=X2),
#                       kappa_Q = kappa_Q, p_Q = p_Q, jumps = jump_counter)
#   })
# }

#' compute G^Q and kappa^Q for Y ~ dist
#' helper function for stressed_sim()
#'
#' @param t A number,time
#' @param x A vector of numbers, state
#' @param eta float
#' @param kappa integer
#' @param stress_type one of "VaR", "VaR & CVaR"
#' @param stress_parms list of parameters for stress
#' VaR: q: stressed VaR value, c: VaR level, VaR_stress: multiplier for P-VaR
#' CVAR: q: stressed VaR value, c: VaR level, s: stressed CVaR level,
#' VaR_stress: multiplier for P-VaR, CVaR_stress: multiplier for P-CVaR
#' @param dist RPS_dist object
#' @param Draws vector
#' @param N_out int, number of draws from distorted jump distribution
#' @param delta_y float, step size for Riemann integration
#' @param y_max int, maximum for integration of h(t,x,y) along y-axes
#' @param tol dbl, cut off integration when prob less than tol
#'
#' @return A list containing:
#' 1) an vector of length(x), containing a draw from G^Q at for each x
#' 2) a vector of length length(x) containing kappa^Q at each x
#'
#' @export
sim_G_kappa_mix <- function(t,x,eta,kappa,stress_type,stress_parms,dist,p,
                        Draws_a,Draws_b,N_out=1,delta_y=0.1,y_max=100,tol=1e-10) {
  with(c(dist,stress_parms), {
    # x is a vector
    xlen <- length(x)

    # find out where to stop riemann sum (when density is below tol)
    y_max_a <- seq(1,y_max)[which(dens_fun_a(seq(1,y_max),parms) < tol)[1]]
    y_max_b <- seq(1,y_max)[which(dens_fun_b(seq(1,y_max),parms) < tol)[1]]
    y_max <- max(y_max_a,y_max_b)
    y <- seq(0,y_max,by=delta_y)
    ylen <- length(y)

    # empty array for h
    h <- array(dim=c(xlen,ylen))

    # compute h at each y
    for (i in 1:ylen){
      h[,i] <- switch(stress_type,
                      "VaR" = h_VaR(x=x,y=y[i],t=t,dist=dist,eta=eta,q=q,kappa=kappa),
                      "CVaR" = h_CVaR(x=x,y=y[i],t=t,dist=dist,eta=eta,q=q,kappa=kappa))
    }

    # initialize
    kappa_Q <- rep(NA,xlen)
    p_Q <- rep(NA,xlen)
    G_Q_1 <- array(dim=c(xlen,N_out))
    G_Q_2 <- array(dim=c(xlen,N_out))

    for (j in 1:xlen){
      # compute integral (kappa.Q):
      h_integral_a <- delta_y * sum(h[j,] * dens_fun_a(y,parms))
      h_integral_b <- delta_y * sum(h[j,] * dens_fun_b(y,parms))
      kappa_Q[j] <- 1 - p + p * h_integral_a
      # new weight between distributions
      p_Q[j] <- p * h_integral_a / kappa_Q[j]
      # draw from G.Q
      u <- rbinom(1,size=1,prob=p_Q[j])
      if (u == 1) weights <- approx(x=y,y=h[j,],xout=Draws_a)$y
      else weights <- approx(x=y,y=h[j,],xout=Draws_b)$y
      weights <- weights / kappa_Q[j]
      G_Q_1[j,] <- u * sample(Draws_a,N_out,replace = TRUE,prob=weights)
      G_Q_2[j,] <- (1-u) * sample(Draws_b,N_out,replace = TRUE,prob=weights)

      # which_dist <- sample(c("a","b"),1,prob=c(p_Q[j],1-p_Q[j]))
      # if (which_dist == "a") {
      #   weights <- approx(x=y,y=h[j,],xout=Draws_a)$y / kappa_Q[j]    # compute weights
      #   G_Q_1[j,] <- sample(Draws_a,N_out,replace = TRUE,prob=weights)   # sample 1 from G^Q_a
      #   G_Q_2[j,] <- 0 # X_2 does not jump!
      # } else if (which_dist == "b") {
      #   weights <- approx(x=y,y=h[j,],xout=Draws_b)$y / kappa_Q[j]    # compute weights
      #   G_Q_1[j,] <- 0 # X_1 does not jump!
      #   G_Q_2[j,] <- sample(Draws_b,N_out,replace = TRUE,prob=weights) # sample 1 from G^Q_b
      # }
    }
    return(list(kappa_Q = kappa * kappa_Q,
                G_Q_1  = G_Q_1,
                G_Q_2 = G_Q_2,
                p_Q = p_Q))
  })
}


sim_baseline_mixture <- function(object=NULL,jump_dist=NULL,
                                 kappa=NULL,Npaths=NULL,time_vec=NULL){
  if (!is.null(object)){
    jump_dist <- object$jump_dist
    Npaths <- dim(object$paths$X1)[2]
    time_vec <- object$time_vec
    kappa <- object$kappa
  }
  # draw from the jump size distribution
  withr::with_seed(720,{
    # N draws from jump size distribution
    Draws_a <- jump_dist$sim_fun_a(1e4,jump_dist$parms)
    Draws_b <- jump_dist$sim_fun_b(1e4,jump_dist$parms)

    Nsteps <- length(time_vec)

    # empty matrices to store results
    X1 <- matrix(nrow=Nsteps,ncol=Npaths)
    X2 <- matrix(nrow=Nsteps,ncol=Npaths)

    X1[1,] <- rep(0,Npaths)
    X2[1,] <- rep(0,Npaths)

    for (i in 2:Nsteps){
      u <- rbinom(Npaths,size=1,prob=jump_dist$parms$p)
      G_1 <- u * sample(Draws_a,Npaths,replace = TRUE)
      G_2 <- (1-u) * sample(Draws_b,Npaths,replace = TRUE)
      U <- runif(Npaths)
      dt <- time_vec[i] - time_vec[i-1]
      X1[i,] <- X1[i-1,] + G_1 * as.integer(U < (1 - exp(-kappa * dt)))
      X2[i,] <- X2[i-1,] + G_2 * as.integer(U < (1 - exp(-kappa * dt)))
      }
    })
    return(list(X1=X1,X2=X2))
}


plot_copula_contour_mixture <- function(object,time){

  time_index <- match(time,object$time_vec)
  u <- seq(0, 1, length.out = 25)
  grid <- as.matrix(expand.grid("x" = u, "y" = u))


  # stressed
  X1 <- object$paths$X1[time_index,]
  X2 <- object$paths$X2[time_index,]
  draws <- cbind(X1,X2)
  ec <-  copula::C.n(grid, X = draws)
  # contourplot2(cbind(as.matrix(grid),ec))
  res_Q <- cbind(as.matrix(grid),ec)
  colnames(res_Q) <- c("x", "y", "z")
  res_Q <- dplyr::as_tibble(cbind(res_Q,measure="stressed")) %>%
    dplyr::mutate(dplyr::across(1:3,as.numeric))
  plot1 <- ggplot2::ggplot(res_Q, ggplot2::aes(x,y,z=z)) +
    ggplot2::geom_contour_filled() +
    ggplot2::theme_minimal()


  # baseline
  baseline <- sim_baseline_mixture(object)
  X1 <- baseline$X1[time_index,]
  X2 <- baseline$X2[time_index,]
  draws <- cbind(X1,X2)
  ec <- copula::C.n(grid, X = draws)
  # contourplot2(cbind(as.matrix(grid),ec))
  res_P <- cbind(as.matrix(grid),ec)
  colnames(res_P) <- c("x", "y", "z")
  res_P <- dplyr::as_tibble(cbind(res_P,measure="original")) %>%
    dplyr::mutate(dplyr::across(1:3,as.numeric))
  plot2 <- ggplot2::ggplot(res_P, ggplot2::aes(x,y,z=z)) +
    ggplot2::geom_contour_filled() +
    ggplot2::theme_minimal()

  xy <- NULL
  for (constant in seq(0.1,0.9,by=0.1)){
    xy <- rbind(xy, dplyr::as_tibble(list(x = seq(0,1,0.01),
                                   y = sapply(seq(0,1,0.01),
                                              function(x) min(constant/x,1)),
                                   constant = constant)))
  }
  xy <-dplyr:: mutate(xy,dplyr::across(3,as.numeric),
               y = ifelse(y==1,NA,y))

  # combined plot
  data <- rbind(res_P,res_Q) %>% dplyr::mutate(measure = forcats::fct_rev(measure))
  plot <- data %>%
    ggplot2::ggplot() +
    ggplot2::geom_contour(ggplot2::aes(x,y,z=z,colour=measure,
                                       linetype = "simulated"),lwd=1) +
    ggplot2::geom_line(data = xy, ggplot2::aes(x,y,group=constant,
                                      linetype = "true independence"), lty = 2) +
    ggplot2::theme_minimal() +
    ggplot2::scale_linetype_manual(values=c("simulated" = 1, "true independence" = 3)) +
    ggplot2::xlab(expression(X[1])) +
    ggplot2::ylab(expression(X[2]))

  return(list(data = data,
              plot = plot))
}


plot_copula_mixture <- function(object,time){
  time_index <- match(time,object$time_vec)

  # stressed
  X1 <- object$paths$X1[time_index,]
  X2 <- object$paths$X2[time_index,]
  stressed <- emp_copula(X1,X2) %>%
    dplyr::mutate(type = "stressed")

  # baseline
  baseline <- sim_baseline_mixture(object)
  X1 <- baseline$X1[time_index,]
  X2 <- baseline$X2[time_index,]
  baseline <- emp_copula(X1,X2) %>%
    dplyr::mutate(type = "baseline")

  browser()

  rbind(baseline,stressed) %>%
    ggplot2::ggplot(ggplot2::aes(U1,U2,colour=type)) +
    ggplot2::facet_grid(~type)+
    ggplot2::geom_point(size=1,alpha=0.4) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::scale_colour_manual(values=c("#00BFC4","#F8766D"))
  # ggplot2::scale_alpha_discrete(range=c(0.7,0.4))
}




plot_copula_dens_contour_mixture <- function(object,time){

  time_index <- match(time,object$time_vec)

  # stressed
  X1 <- object$paths$X1[time_index,]
  X2 <- object$paths$X2[time_index,]
  u1 <- copula::pobs(cbind(X1,X2))
  kde.fit1 <- kdecopula::kdecop(u1)
  plot(kde.fit1)
  graphics::contour(kde.fit1,margins="unif",#levels = c(0.1,0.5,1.5,2,3),
                    col="#F8766D",lwd=2,cex.lab=1.5,cex.axis=1.5,
                    labcex=1.1,vfont=c("sans serif", "bold"))

  # baseline
  baseline <- sim_baseline_mixture(object)
  X1 <- baseline$X1[time_index,]
  X2 <- baseline$X2[time_index,]
  u2 <- copula::pobs(cbind(X1,X2))
  kde.fit2 <- kdecopula::kdecop(u2)
  plot(kde.fit2)
  graphics::contour(kde.fit2,margins="unif",#levels = c(0.1,0.5,1.5,2,3),
                    col="#00BFC4",add=TRUE,lwd=2,
                    labcex=1.1,vfont=c("sans serif", "bold"))
  # add legend
  graphics::legend("bottomright", legend = c("Q", "P"),
                   col = c("#F8766D","#00BFC4"),lty=1:1,lwd=c(2,2),cex=1.2)

}

#
# baseline_ex <- sim_baseline_mixture(mixed_ex3)
# ks.test(baseline_ex$X2[501,],mixed_ex3$paths$X2[501,])
# pvec <- rep(NA,length(mixed_ex3$time_vec) )
# for (i in 1:length(mixed_ex3$time_vec)){
#   pvec[i] <- ks.test(baseline_ex$X2[i,],mixed_ex3$paths$X2[i,])$p.value
# }
# plot(mixed_ex2$time_vec,pvec)
#
# time <- mixed_ex3$time_vec
# as_tibble(cbind(time,pvec)) %>% ggplot(aes(time,pvec)) + geom_point() +
#   geom_hline(yintercept = 0.05)
#
# baseline_ex <- sim_baseline_mixture(mixed_ex_b)
# ks.test(baseline_ex$X2[101,],mixed_ex_b$paths$X2[101,])
#
# ks.test(rnorm(1e4),rnorm(1e4))
# ks.test(rnorm(1e4),rgamma(1e4,2,1))

#===============================================================================
#
# # density function
# dMixture <- function(x,p,alpha,beta){
#   p * dgamma(x,shape=alpha,rate=beta) + (1-p) * ifelse(x==0, 1, 0)
# }
#
# # distribution function
# pMixture <- function(q,p,alpha,beta){
#   p * pgamma(q,shape=alpha,rate=beta) + (1-p) * ifelse(q >= 0, 1, 0)
# }
#
# # simulation function
# rMixture <- function(n,p,alpha,beta){
#   u <- rbinom(n,size=1,prob=p)
#   u *  rgamma(n,shape=alpha,rate=beta) + (1-u) * 0
# }
#
# # characteristic function
# cMixture <- function(x,p,alpha,beta){
#   p * (1 - 1i*x/beta)^(-alpha) + 1 - p
# }

#------
#
# gamma_mixed <- new_RPS_dist_mix(sim_fun_a = function(x,parms) rgamma(x,shape=parms$alpha1,rate=parms$beta1),
#   sim_fun_b = function(x,parms) rgamma(x,shape=parms$alpha2,rate=parms$beta2),
#   dens_fun_a = function(x,parms) dgamma(x,shape=parms$alpha1,rate=parms$beta1),
#   dens_fun_b = function(x,parms) dgamma(x,shape=parms$alpha2,rate=parms$beta2),
#   char_fun = function(x,parms) parms$p * (1 - 1i*x/parms$beta1)^(-parms$alpha1) + 1 - parms$p,
#   mean_fun = function(parms) parms$p*parms$alpha1/parms$beta1,
#   parms = list(alpha1 = 2,beta1=1,alpha2=2,beta2=1,p = 0.5))
# #
# tm <- proc.time()
# mixed_ex <- stressed_sim_mix(kappa = 5, jump_dist = gamma_mixed,
#                              stress_type = "VaR",
#                              stress_parms = list(c=0.9,VaR_stress=1.15),
#                              Npaths=1e4, endtime=1, dt=2e-2)
# proc.time() - tm
#
#
# tm <- proc.time()
# mixed_ex3 <- stressed_sim_mix(kappa = 5, jump_dist = gamma_mixed,
#                               stress_type = "VaR",
#                               stress_parms = list(c=0.9,VaR_stress=1.15),
#                               Npaths=1e4, endtime=1, dt=2e-3)
# proc.time() - tm
#
# tm <- proc.time()
# mixed_ex <- stressed_sim_mix(kappa = 5, jump_dist = gamma_mixed,
#                                   stress_type = "VaR",
#                                   stress_parms = list(c=0.9,VaR_stress=1.1),
#                                   Npaths=1e4, endtime=1, dt=1e-3)
# proc.time() - tm


