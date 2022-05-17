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
#' @param end_time float, when to end sim
#' @param dt float, step size in time
#'
#' @return RPS_model object
#' @export
stressed_sim_biv <- function(kappa, jump_dist, stress_type = "VaR",
                         stress_parms = list(c=c,q=q,s=s,VaR_stress=VaR_stress,
                                             CVaR_stress=CVaR_stress),
                         Npaths=1e4, end_time=1, dt=1e-2, beep = TRUE){

  # first, draw from the jump size distribution
  Ndraws <- 1e4
  with(jump_dist,{
    withr::with_seed(720,{

      Draws <- copula::rMvdc(Ndraws,biv_dist)  # N draws from jump size distribution

      Nsteps <- end_time / dt + 1   # number of steps to take

      # if time_stress is empty, assume stress is at the end
      if (is.null(stress_parms$time_stress)) {
        cat("No stress time specified. Applying stress at terminal time. \n")
        stress_parms$time_stress <- end_time
      }

      if (is.null(stress_parms$q) & !(is.null(stress_parms$VaR_stress))){
        stress_parms$q <- stress_parms$VaR_stress *
          compute_VaR(stress_parms$time_stress,kappa,stress_parms$c,jump_dist)
      }
      if ((stress_type == "CVaR") & is.null(stress_parms$s) &
          !(is.null(stress_parms$CVaR_stress))){
        stress_parms$s <- stress_parms$CVaR_stress *
          compute_CVaR(stress_parms$time_stress,kappa=kappa,
                       q=stress_parms$q,c=stress_parms$c,dist=jump_dist)
        if (stress_parms$s < stress_parms$q) {
          stop("Incompatible parameter choice: attempt to set CVaR below VaR.\n")
        }
      }

      eta <- switch(stress_type,
                    "VaR" = eta_VaR(kappa=kappa, stress_parms=stress_parms,
                                    dist=jump_dist),
                    "CVaR" = eta_CVaR(kappa=kappa, stress_parms=stress_parms,
                                      dist=jump_dist))
      cat('eta:',eta, "\n")

      times <- seq(0,end_time,by=dt)

      # initialize
      X1 <- matrix(nrow=Nsteps,ncol=Npaths) # empty matrix to store results
      X2 <- matrix(nrow=Nsteps,ncol=Npaths)
      kappa_Q <- matrix(nrow=Nsteps,ncol=Npaths)

      X1[1,] <- rep(0,Npaths)   # X starts at 0 for all paths
      X2[1,] <- rep(0,Npaths)
      times <- seq(0,end_time,by=dt)

      for (i in 2:Nsteps){
        U <- stats::runif(Npaths)  # generate Npaths independent unif[0,1]

        # compute directly
        distorted <- sim_G_kappa_biv(t=times[i], x=X1[i-1,], eta=eta,
                                     kappa=kappa, stress_type=stress_type,
                                     stress_parms=stress_parms,
                                     dist=jump_dist, Draws=Draws)
        kappa_Q[i,] <- distorted$kappa_Q  # store kappa
        G_Q_draw1 <- distorted$G_Q_1
        G_Q_draw2 <- distorted$G_Q_2

        # simulate forward:
        jump <- as.integer(U < (1 - exp(-kappa_Q[i,] * dt)))
        X1[i,] <- X1[i-1,] + G_Q_draw1 * jump
        X2[i,] <- X2[i-1,] + G_Q_draw2 * jump
      }
    })
    if (beep == TRUE) beepr::beep()
    new_RPS_model(model_type = "bivariate_copula", jump_dist = jump_dist,
                  kappa = kappa, stress_type = stress_type,
                  stress_parms = stress_parms, time_vec = times,
                  eta = eta, paths = list(X1=X1, X2=X2), intensity = kappa_Q)
  })
}

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
sim_G_kappa_biv <- function(t,x,eta,kappa,stress_type,stress_parms,dist,
                        Draws,N_out=1,delta_y=0.1,y_max=100,tol=1e-10) {
  with(c(dist,stress_parms), {
    # x is a vector
    xlen <- length(x)
    y_max <- seq(1,y_max)[which(dens_fun(seq(1,y_max),parms) < tol)[1]]
    y <- seq(0,y_max,by=delta_y)
    ylen <- length(y)
    h <- array(dim=c(xlen,ylen))

    # compute h at each y
    for (i in 1:ylen){
      h[,i] <- switch(stress_type,
                      "VaR" = h_VaR(x=x,y=y[i],t=t,dist=dist,eta=eta,
                                    kappa=kappa,stress_parms=stress_parms),
                      "CVaR" = h_CVaR(x=x,y=y[i],t=t,dist=dist,eta=eta,
                                      kappa=kappa,stress_parms=stress_parms))
    }

    # initialize
    kappa_Q <- rep(NA,xlen)
    G_Q_1 <- array(dim=c(xlen,N_out))
    G_Q_2 <- array(dim=c(xlen,N_out))

    for (j in 1:xlen){
      # compute integral (kappa.Q):
      kappa_Q[j] <- delta_y * sum(h[j,] * dens_fun(y,parms))
      # draw from G.Q
      weights <- approx(x=y,y=h[j,],xout=Draws[,1])$y / kappa_Q[j]    # compute weights
      G_Q_1[j,] <- sample(Draws[,1],N_out,replace = TRUE,prob=weights)   # sample 1 from G^Q_1
      G_Q_2[j,] <- sample(Draws[,2],N_out,replace = TRUE,prob=weights)   # sample 1 from G^Q_2
    }
    return(list(kappa_Q = kappa * kappa_Q,
                G_Q_1  = G_Q_1,
                G_Q_2 = G_Q_2))
  })
}


#' function to get path of X under original probability measure
#'
#' @param object RPS_model object
#'
#' @return matrix, paths of original model
#' @export
#'
sim_baseline <- function(object){
  with(object,{
    # draw from the jump size distribution
    withr::with_seed(720,{
      Draws <- jump_dist$sim_fun(1e4,jump_dist$parms)

      Npaths <- dim(paths)[2]
      Nsteps <- length(time_vec)

      X <- matrix(nrow=Nsteps,ncol=Npaths) # empty matrix to store results
      X[1,] <- rep(0,Npaths)

      for (i in 2:Nsteps){
        U <- stats::runif(Npaths)
        dt <- time_vec[i] - time_vec[i-1]
        X[i,] <- X[i-1,] + sample(Draws,Npaths,replace=TRUE) * as.integer(U < (1 - exp(-kappa * dt)))
      }
    })
    return(X)
  })
}

sim_baseline_biv <- function(object=NULL,jump_dist=NULL,
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
      Draws <- copula::rMvdc(1e4,jump_dist$biv_dist)

      Nsteps <- length(time_vec)

      # empty matrices to store results
      X1 <- matrix(nrow=Nsteps,ncol=Npaths)
      X2 <- matrix(nrow=Nsteps,ncol=Npaths)

      X1[1,] <- rep(0,Npaths)
      X2[1,] <- rep(0,Npaths)

      for (i in 2:Nsteps){
        U <- stats::runif(Npaths)
        dt <- time_vec[i] - time_vec[i-1]
        X1[i,] <- X1[i-1,] + sample(Draws[,1],Npaths,replace=TRUE) *
          as.integer(U < (1 - exp(-kappa * dt)))
        X2[i,] <- X2[i-1,] + sample(Draws[,2],Npaths,replace=TRUE) *
          as.integer(U < (1 - exp(-kappa * dt)))
      }
    })
    return(list(X1=X1,X2=X2))
}



### MIXTURE

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
#' @param end_time float, when to end sim
#' @param dt float, step size in time
#'
#' @return RPS_model object
#' @export
stressed_sim_mix <- function(kappa, jump_dist, stress_type = "VaR",
                             stress_parms = list(c=c,q=q,s=s,
                                                 VaR_stress=VaR_stress,
                                                 CVaR_stress=CVaR_stress),
                             Npaths=1e4,end_time=1, dt=1e-2){

  # first, draw from the jump size distribution
  Ndraws <- 1e4
  with(jump_dist,{
    withr::with_seed(720,{
      # N draws from both jump size distributions
      Draws_a <- sim_fun_a(Ndraws,parms)
      Draws_b <- sim_fun_b(Ndraws,parms)

      Nsteps <- end_time / dt + 1   # number of steps to take
      times <- seq(0,end_time,by=dt)

      # if time_stress is empty, assume stress is at the end
      if (is.null(stress_parms$time_stress)) {
        cat("No stress time specified. Applying stress at terminal time. \n")
        stress_parms$time_stress <- end_time
      }

      if (is.null(stress_parms$q) & !(is.null(stress_parms$VaR_stress))){
        stress_parms$q <- stress_parms$VaR_stress *
          compute_VaR(stress_parms$time_stress,kappa,stress_parms$c,jump_dist)
      }
      if ((stress_type == "CVaR") & is.null(stress_parms$s) &
          !(is.null(stress_parms$CVaR_stress))){
        stress_parms$s <- stress_parms$CVaR_stress *
          compute_CVaR(stress_parms$time_stress,kappa=kappa,
                       q=stress_parms$q,c=stress_parms$c,dist=jump_dist)
        if (stress_parms$s < stress_parms$q) {
          stop("Incompatible parameter choice: attempt to set CVaR below VaR.\n")
        }
      }

      eta <- switch(stress_type,
                    "VaR" = eta_VaR(kappa=kappa, stress_parms=stress_parms,
                                    dist=jump_dist),
                    "CVaR" = eta_CVaR(kappa=kappa, stress_parms=stress_parms,
                                      dist=jump_dist))

      # initialize
      X1 <- matrix(nrow=Nsteps,ncol=Npaths) # empty matrix to store results
      X2 <- matrix(nrow=Nsteps,ncol=Npaths)
      kappa_Q <- matrix(nrow=Nsteps,ncol=Npaths)
      p_Q <- matrix(nrow=Nsteps,ncol=Npaths)

      X1[1,] <- rep(0,Npaths)   # X starts at 0 for all paths
      X2[1,] <- rep(0,Npaths)
      times <- seq(0,end_time,by=dt)

      for (i in 2:Nsteps){
        distorted <- sim_G_kappa_mix(t=times[i],x=X1[i-1,],eta=eta,
                                     kappa=kappa,stress_type=stress_type,
                                     stress_parms=stress_parms,p=parms$p,
                                     dist=jump_dist,Draws_a=Draws_a,
                                     Draws_b=Draws_b)
        U <- stats::runif(Npaths)  # generate Npaths independent unif[0,1]

        # compute kappa.Q and draw from G.Q, using interpolation
        kappa_Q[i,] <- distorted$kappa_Q  # store kappa
        p_Q[i,] <- distorted$p_Q # store p
        G_Q_draw1 <- distorted$G_Q_1
        G_Q_draw2 <- distorted$G_Q_2

        # simulate forward:
        jump <- as.integer(U < (1 - exp(-kappa_Q[i,] * dt)))
        X1[i,] <- X1[i-1,] + G_Q_draw1 * jump
        X2[i,] <- X2[i-1,] + G_Q_draw2 * jump

      }
    })
    beepr::beep()
    new_RPS_model(model_type = "bivariate_mixture", jump_dist = jump_dist,
                  kappa = kappa, stress_type = stress_type,
                  stress_parms = stress_parms, time_vec = times,
                  eta = eta, paths = list(X1=X1, X2=X2),
                  intensity = list(kappa_Q=kappa_Q, p_Q=p_Q))
  })
}



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
                      "VaR" = h_VaR(x=x,y=y[i],t=t,dist=dist,eta=eta,
                                    kappa=kappa,stress_parms=stress_parms),
                      "CVaR" = h_CVaR(x=x,y=y[i],t=t,dist=dist,eta=eta,
                                      kappa=kappa,stress_parms=stress_parms))
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
      u <- stats::rbinom(1,size=1,prob=p_Q[j])
      if (u == 1) weights <- approx(x=y,y=h[j,],xout=Draws_a)$y
      else weights <- approx(x=y,y=h[j,],xout=Draws_b)$y
      weights <- weights / kappa_Q[j]
      G_Q_1[j,] <- u * sample(Draws_a,N_out,replace = TRUE,prob=weights)
      G_Q_2[j,] <- (1-u) * sample(Draws_b,N_out,replace = TRUE,prob=weights)

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
      u <- stats::rbinom(Npaths,size=1,prob=jump_dist$parms$p)
      G_1 <- u * sample(Draws_a,Npaths,replace = TRUE)
      G_2 <- (1-u) * sample(Draws_b,Npaths,replace = TRUE)
      U <- stats::runif(Npaths)
      dt <- time_vec[i] - time_vec[i-1]
      X1[i,] <- X1[i-1,] + G_1 * as.integer(U < (1 - exp(-kappa * dt)))
      X2[i,] <- X2[i-1,] + G_2 * as.integer(U < (1 - exp(-kappa * dt)))
    }
  })
  return(list(X1=X1,X2=X2))
}
