#' Simulate the stressed model
#'
#' @param kappa double, compound Poisson parameter
#' @param jump_dist RPS_dist object, contains jump size distribution
#' @param stress_type one of "VaR", "VaR & CVaR"
#' @param stress_parms list of parameters for stress
#' VaR: q: stressed VaR value, c: VaR level, VaR_stress: multiplier for P-VaR
#' CVAR: q: stressed VaR value, c: VaR level, s: stressed CVaR level,
#' VaR_stress: multiplier for P-VaR, CVaR_stress: multiplier for P-CVaR
#' time_stress: time at which stress occurs; if NULL, assumed to be end_time
#' @param Npaths integer, number of paths
#' @param end_time float, when to end sim
#' @param dt float, step size in time
#'
#' @return RPS_model object
#' @export
stressed_sim <- function(kappa, jump_dist, stress_type = "VaR",
                         stress_parms = list(c=c,q=q,s=s,
                                             VaR_stress=VaR_stress,
                                             CVaR_stress=CVaR_stress,
                                             time_stress = time_stress),
                         Npaths=1e4, end_time=1, dt=1e-2){

  # check type of jump distribution
  stopifnot(is(jump_dist) == "RPS_dist_univ")
  # check parameters
  stopifnot(dt < end_time)
  stopifnot(stress_parms$time_stress <= end_time)

  # first, draw from the jump size distribution
  Ndraws <- 1e4
  with(jump_dist,{
    withr::with_seed(720,{
      Draws <- sim_fun(Ndraws,parms)  # N draws from jump size distribution

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

      cat("eta:", eta, "\n")


      # Initializes the progress bar
      pb <- txtProgressBar(min = 0, max = Nsteps, style = 3,
                           width = 50, char = "=")

      # create grid of Gs and kappas
      times <- seq(0,end_time,by=dt)

      # initialize
      X <- matrix(nrow=Nsteps,ncol=Npaths) # empty matrix to store results
      kappa_Q <- matrix(nrow=Nsteps,ncol=Npaths)

      X[1,] <- rep(0,Npaths)   # X starts at 0 for all paths
      times <- seq(0,end_time,by=dt)

      for (i in 2:Nsteps){
        U <- runif(Npaths)  # generate Npaths independent unif[0,1]

        # compute kappa.Q and draw from G.Q
        distorted <- sim_G_kappa(t=times[i],x=X[i-1,],eta=eta,kappa=kappa,
                                 stress_type=stress_type,
                                 stress_parms=stress_parms,
                                 dist=jump_dist,Draws=Draws)
        kappa_Q[i,] <- distorted$kappa_Q
        G_Q_draw <- distorted$G_Q

        # simulate forward:
        X[i,] <- X[i-1,] + G_Q_draw * as.integer(U < (1 - exp(-kappa_Q[i,] * dt)))

        # Sets the progress bar to the current state
        setTxtProgressBar(pb, i)
        }
      })
    beepr::beep()
    RPS_model(model_type = "univariate", jump_dist = jump_dist,
              kappa = kappa, stress_type = stress_type,
              stress_parms = stress_parms, time_vec = times,
              eta = eta, paths = X, intensity = kappa_Q)
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
#'
#' @return A list containing:
#' 1) an vector of length(x), containing a draw from G^Q at for each x
#' 2) a vector of length length(x) containing kappa^Q at each x
#'
#' @export
sim_G_kappa <- function(t,x,eta,kappa,stress_type,stress_parms,dist,
                        Draws,N_out=1,delta_y=0.1) {
  with(c(dist,stress_parms), {
    # x is a vector
    xlen <- length(x)
    y <- seq(dist$min,dist$max,by=delta_y)
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
    G_Q <- array(dim=c(xlen,N_out))

    for (j in 1:xlen){
      # compute integral (kappa.Q):
      kappa_Q[j] <- delta_y * sum(h[j,] * dens_fun(y,parms))
      # draw from G.Q
      weights <- approx(x=y,y=h[j,],xout=Draws)$y / kappa_Q[j]    # compute weights
      weights <- ifelse(weights < 0,0,weights)
      G_Q[j,] <- sample(Draws,N_out,replace = TRUE,prob=weights)   # sample 1 from G^Q
    }
    return(list(kappa_Q = kappa * kappa_Q,
                G_Q = G_Q))
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
        U <- runif(Npaths)
        dt <- time_vec[i] - time_vec[i-1]
        X[i,] <- X[i-1,] + sample(Draws,Npaths,replace=TRUE) * as.integer(U < (1 - exp(-kappa * dt)))
      }
    })
    return(X)
  })
}
