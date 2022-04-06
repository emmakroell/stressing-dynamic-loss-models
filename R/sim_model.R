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
stressed_sim <- function(kappa, jump_dist, stress_type = "VaR",
                         stress_parms = list(c=c,q=q,s=s,VaR_stress=VaR_stress,
                                             CVaR_stress=CVaR_stress),
                         Npaths=1e4,endtime=1, dt=1e-2){

  # first, draw from the jump size distribution
  Ndraws <- 1e4
  with(jump_dist,{
    withr::with_seed(720,{
      Draws <- copula::rMvdc(Ndraws,biv_dist)  # N draws from jump size distribution

      Nsteps <- endtime / dt + 1   # number of steps to take

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
      cat('eta:',eta)

      # create grid of Gs and kappas
      times <- seq(0,endtime,by=dt)
      # choose X max based on parameters:
      x_max <- 8 * kappa * mean(jump_dist) # FIND A BETTER WAY TO DO THIS
      x1_seq <- seq(0,x_max,by=2e-2)
      x2_seq <- seq(0,x_max,by=2e-2)
      len_x <- length(x1_seq)

      G_grid1 <- array(dim=c(Nsteps,len_x))
      G_grid2 <- array(dim=c(Nsteps,len_x))
      kappa_grid <- matrix(nrow = Nsteps, ncol = len_x)


      # compute G and kappa at grid points
      for (i in 1:Nsteps){
        distorted <- sim_G_kappa(t=times[i],x=x1_seq,eta=eta,kappa=kappa,
                                 stress_type=stress_type,
                                 stress_parms=stress_parms,
                                 dist=jump_dist,Draws=Draws)

        G_grid1[i,] <- distorted$G_Q_1
        G_grid2[i,] <- distorted$G_Q_2
        kappa_grid[i,] <- distorted$kappa_Q
      }

      # initialize
      X1 <- matrix(nrow=Nsteps,ncol=Npaths) # empty matrix to store results
      X2 <- matrix(nrow=Nsteps,ncol=Npaths)
      kappa_Q <- matrix(nrow=Nsteps,ncol=Npaths)

      X1[1,] <- rep(0,Npaths)   # X starts at 0 for all paths
      X2[1,] <- rep(0,Npaths)
      times <- seq(0,endtime,by=dt)

      for (i in 2:Nsteps){
        U <- runif(Npaths)  # generate Npaths independent unif[0,1]

        # compute kappa.Q and draw from G.Q, using interpolation
        kappa_Q[i,] <- approx(x=x1_seq,y=kappa_grid[i-1,],xout=X1[i-1,])$y  # store kappa
        G_Q_draw1 <- approx(x=x1_seq,y=G_grid1[i-1,],xout=X1[i-1,])$y
        G_Q_draw2 <- approx(x=x2_seq,y=G_grid2[i-1,],xout=X2[i-1,])$y
        # simulate forward:
        X1[i,] <- X1[i-1,] + G_Q_draw1 * as.integer(U < (1 - exp(-kappa_Q[i,] * dt)))
        X2[i,] <- X2[i-1,] + G_Q_draw2 * as.integer(U < (1 - exp(-kappa_Q[i,] * dt)))

        # if X passed the grid max, return an error:
        if (any(c(X1[i,],X2[i,]) > x_max)){
          cat("max X1 ", max(X1[i,]), "max X2 ", max(X2[i,]), "\n")
          stop("Path value exceeded grid max. Increase x max. \n")
          }
      }
    })
    beepr::beep()
    new_RPS_model(jump_dist = jump_dist, kappa = kappa,
                  stress_type = stress_type, stress_parms = stress_parms,
                  time_vec = times, paths = list(X1=X1, X2=X2),
                  kappa_Q = kappa_Q)
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
sim_G_kappa <- function(t,x,eta,kappa,stress_type,stress_parms,dist,
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
                      "VaR" = h_VaR(x=x,y=y[i],t=t,dist=dist,eta=eta,q=q,kappa=kappa),
                      "CVaR" = h_CVaR(x=x,y=y[i],t=t,dist=dist,eta=eta,q=q,kappa=kappa))
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
      # weights <- ifelse(weights < 0,0,weights)
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
        U <- runif(Npaths)
        dt <- time_vec[i] - time_vec[i-1]
        X[i,] <- X[i-1,] + sample(Draws,Npaths,replace=TRUE) * as.integer(U < (1 - exp(-kappa * dt)))
      }
    })
    return(X)
  })
}

sim_baseline_biv <- function(object){
  with(object,{
    # draw from the jump size distribution
    withr::with_seed(720,{
      # N draws from jump size distribution
      Draws <- copula::rMvdc(1e4,jump_dist$biv_dist)

      Npaths <- dim(paths$X1)[2]
      Nsteps <- length(time_vec)

      # empty matrices to store results
      X1 <- matrix(nrow=Nsteps,ncol=Npaths)
      X2 <- matrix(nrow=Nsteps,ncol=Npaths)

      X1[1,] <- rep(0,Npaths)
      X2[1,] <- rep(0,Npaths)

      for (i in 2:Nsteps){
        U <- runif(Npaths)
        dt <- time_vec[i] - time_vec[i-1]
        X1[i,] <- X1[i-1,] + sample(Draws[,1],Npaths,replace=TRUE) *
          as.integer(U < (1 - exp(-kappa * dt)))
        X2[i,] <- X2[i-1,] + sample(Draws[,2],Npaths,replace=TRUE) *
          as.integer(U < (1 - exp(-kappa * dt)))
      }
    })
    return(list(X1=X1,X2=X2))
  })
}
