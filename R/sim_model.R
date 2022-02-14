#' Simulate the stressed model
#'
#' @param kappa double, compound Poisson parameter
#' @param jump_dist RPS_dist object, contains jump size distribution
#' @param stress_type one of "VaR", "VaR & CVaR"
#' @param stress_parms list of parameters for stress
#' VaR: q: stressed VaR value, c: VaR level
#' CVAR: TBA
#' @param Npaths integer, number of paths
#' @param endtime float, when to end sim
#' @param dt float, step size in time
#'
#' @return RPS_model object
#' @export
stressed_sim <- function(kappa, jump_dist, stress_type = "VaR",
                         stress_parms = list(c=c,q=q), Npaths=1e4,
                         endtime=1, dt=1e-2){

  # first, draw from the jump size distribution
  Ndraws <- 1e4
  with(jump_dist,{
    ptm <- proc.time()
    withr::with_seed(720,{
      Draws <- sim_fun(Ndraws,parms)  # N draws from jump size distribution
    })

    Nsteps <- endtime / dt + 1   # number of steps to take

    eta <- VaR_eta(kappa=5, stress_parms=stress_parms, dist=gamma_2_1)

    # create grid of Gs and kappas
    times <- seq(0,endtime,by=dt)
    # choose X max based on parameters:
    x_max <- 5 * kappa * mean(jump_dist)
    x_seq <- seq(0,x_max,by=0.1)
    len_x <- length(x_seq)

    G_grid <- matrix(nrow = Nsteps, ncol = len_x)
    kappa_grid <- matrix(nrow = Nsteps, ncol = len_x)


    # compute G and kappa at grid points
    for (i in 1:Nsteps){
      distorted <- sim_G_kappa(t=times[i],x=x_seq,eta=eta,kappa=kappa,
                               q=stress_parms$q,c=stress_parms$c,
                               dist=jump_dist,Draws=Draws,N=Ndraws)
      G_grid[i,] <- distorted$G_Q
      kappa_grid[i,] <- distorted$kappa_Q
    }

    # initialize
    X <- matrix(nrow=Nsteps,ncol=Npaths) # empty matrix to store results
    kappa_Q <- matrix(nrow=Nsteps,ncol=Npaths)

    X[1,] <- rep(0,Npaths)   # X starts at 0 for all paths
    times <- seq(0,endtime,by=dt)

    for (i in 2:Nsteps){
      U <- runif(Npaths)  # generate Npaths independent unif[0,1]

      # compute kappa.Q and draw from G.Q, using interpolation
      kappa_Q[i,] <- approx(x=x_seq,y=kappa_grid[i-1,],xout=X[i-1,])$y  # store kappa
      G_Q_draw <- approx(x=x_seq,y=G_grid[i-1,],xout=X[i-1,])$y
      # simulate forward:
      X[i,] <- X[i-1,] + G_Q_draw * as.integer(U < (1 - exp(-kappa_Q[i,] * dt)))

      # if X passed the grid max, return an error:
      if (any(X[i,] > x_max)) stop("Path value exceeded grid max. Increase x.max. \n")
    }

    print(proc.time() - ptm)

    new_RPS_model(jump_dist = jump_dist, kappa = kappa,
                  stress_type = stress_type, stress_parms = stress_parms,
                  time_vec = times, paths = X, kappa_Q = kappa_Q)

  })
}


#' compute G^Q and kappa^Q for Y ~ dist
#' helper function for stressed_sim()
#'
#' @param t A number,time
#' @param x A vector of numbers, state
#' @param eta float
#' @param kappa integer
#' @param q float, distorted quantile
#' @param c float, qantile to distort, c in (0,1)
#' @param dist RPS_dist object
#' @param Draws vector
#' @param delta_y float, step size for Riemann integration
#' @param N integer, number of draws from distribution
#'
#' @return A list containing:
#' 1) an vector of length(x), containing a draw from G^Q at for each x
#' 2) a vector of length length(x) containing kappa^Q at each x
#'
#' @export
sim_G_kappa <- function(t,x,eta,kappa,q,c,dist,Draws,delta_y=0.1,N=1e4){
  with(dist, {
    # x is a vector
    xlen <- length(x)
    y <- seq(0,q+5,by=delta_y)  # h is constant for x+y > q
    ylen <- length(y)
    h <- array(dim=c(xlen,ylen))

    # compute h at each y
    for (i in 1:ylen){
      h[,i] <- h_FST(x=x,y=y[i],t=t,dist=dist,eta=eta,q=q,c=c,kappa=kappa)
    }

    # initialize
    kappa_Q <- rep(NA,xlen)
    G_Q <- rep(NA,xlen)

    for (j in 1:xlen){
      # compute integral (kappa.Q):
      kappa_Q[j] <- dist_fun(0,parms) + delta_y * sum(h[j,] * dens_fun(y,parms)) +
        tail(h[j,],1) * (1-dist_fun(tail(y,1),parms))
      # draw from G.Q
      weights <- approx(x=y,y=h[j,],xout=Draws)$y / kappa_Q[j]    # compute weights
      G_Q[j] <- sample(Draws, 1, replace = TRUE, prob=weights)   # sample 1 from G^Q
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
    })

    Npaths <- dim(paths)[2]

    Nsteps <- length(time_vec)
    X <- matrix(nrow=Nsteps,ncol=Npaths) # empty matrix to store results
    X[1,] <- rep(0,Npaths)

    for (i in 2:Nsteps){
      U <- runif(Npaths)
      dt <- time_vec[i] - time_vec[i-1]
      X[i,] <- X[i-1,] + sample(Draws,Npaths) * as.integer(U < (1 - exp(-kappa * dt)))
    }
    return(X)
  })
}
