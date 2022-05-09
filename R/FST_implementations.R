# Functions that use FST

#-------------------------------------------------------------------------------
# VaR
#-------------------------------------------------------------------------------

#' Computes the P probability that X is less than q using FST
#'
#' @param x dbl
#' @param t dbl
#' @param kappa dbl, compound Poisson parameter
#' @param q dbl
#' @param dist RPS_dist object
#' @param grid_max int
#' @param time_stress dbl
#' @param N int
#'
#' @return
#' @export
#'
P_prob <- function(x,t,kappa,q,dist,grid_max=50,time_stress=1,N=8192){
  indic_func <- function(x) as.integer(x < q)
  FST(x=x, t=t, dist=dist, kappa = kappa,terminal_cond=indic_func,
      grid_max=grid_max, time_stress=time_stress, N=N)
}


#' computes VaR under the P measure using FST and stats::uniroot
#'
#' @param time dbl, time at which to compute VaR
#' @param kappa dbl, compound Poisson parameter
#' @param c dbl, level for VaR between 0 and 1
#' @param dist RPS_dist object
#'
#' @return
#' @export
#'
compute_VaR <- function(time,kappa,c,dist){
  grid_max <- 5 * kappa * mean(dist)
  stats::uniroot(function(q) P_prob(x=0,t=0,kappa=kappa,q=q,
                                    dist=dist,grid_max=grid_max,
                                    time_stress=time) - c,
                 interval=c(0, 5*kappa*mean(dist)))$root
}


#' Computes the optimal eta for VaR stress using FST
#'
#' @param kappa double, compound Poisson parameter
#' @param stress_parms list of parameters for stress; must contain:
#' q: stressed VaR value, time_stress: time at which stress occurs
#' @param dist RPS_dist object, contains jump size distribution
#'
#' @return double, optimal eta for VaR stress
#' @export
#'
eta_VaR <- function(kappa,stress_parms,dist) {
  with(stress_parms,{
    grid_max <- 5 * kappa * mean(dist)
    p <- P_prob(x=0, t=0, kappa = kappa, q=q, dist=dist,
                grid_max=grid_max, time_stress=time_stress)
    eta <- log((1-c)/c * p/(1-p))
    return(eta)
  })
}


#' Computes h(t,x,y) for VaR constraint
#'
#' @param x vector
#' @param y vector
#' @param t double
#' @param dist RPS_dist object, contains jump size distribution
#' @param eta double, optimal eta for VaR stress
#' @param q double, stressed VaR value
#' @param kappa double, compound Poisson parameter
#' @param stress_parms list of parameters for stress; must contain:
#' q: stressed VaR value, time_stress: time at which stress occurs
#' @param N int, number of points for FST method
#'
#' @return
#' @export
#'
#'
h_VaR <- function(x, y, t, dist, eta, kappa, stress_parms, N=8192) {
  with(stress_parms,{
    if (t > time_stress) h <- 1
    else{
      grid_max <- max(x) + max(y) + 50
      num <- 1 + (exp(-eta) - 1) * P_prob(x=x+y,t=t,kappa=kappa,q=q,
                                          dist=dist,grid_max=grid_max,
                                          time_stress=time_stress,N=N)
      denom <- 1 + (exp(-eta) - 1) * P_prob(x=x,t=t,kappa=kappa,q=q,
                                            dist=dist,grid_max=grid_max,
                                            time_stress=time_stress,N=N)
      h <- num / denom
    }
    return(h)
  })
}


#-------------------------------------------------------------------------------
# CVaR
#-------------------------------------------------------------------------------


#' Compute CVaR under the P measure
#'
#' @param time dbl, time at which to compute CVaR
#' @param kappa dbl, compound Poisson parameter
#' @param q dbl, stressed VaR value
#' @param c dbl, VaR level, between 0 and 1
#' @param dist RPS_dist object
#'
#' @return
#' @export
#'
compute_CVaR <- function(time,kappa,q,c,dist) {
  terminal_eqn <- function(x) (x - q) * as.integer(x > q)
  grid_max <- 5 * kappa * mean(dist)
  expectation <- FST(x=0, t=0, dist=dist, kappa = kappa,
                     terminal_cond=terminal_eqn,
                     grid_max=grid_max, time_stress=time)
  1/(1-c) * expectation + q
}


#' Computes P-expectation of e^(-eta2 * (x-q)) indicator(x >= q)
#' helper function for eta_CVaR and h_CVaR
#'
#' @param x dbl
#' @param t dbl
#' @param kappa dbl, compound Poisson parameter
#' @param q dbl, stressed VaR value
#' @param eta2 dbl
#' @param dist RPS_dist object
#' @param grid_max int
#' @param time_stress dbl
#' @param N int
#'
#' @return
#' @export
#'
P_expectation <- function(x,t,kappa,q,eta2,dist,grid_max,time_stress=1,N=8192){
  terminal_eqn <- function(x) exp(-eta2 * (x - q)) * as.integer(x >= q)
  FST(x=x, t=t, dist=dist, kappa=kappa, terminal_cond=terminal_eqn,
      grid_max=grid_max, time_stress=time_stress, N=N)
}


#' Equation to be solved to find eta_2
#'
#' @param eta dbl
#' @param dist RPS_dist object
#' @param q dbl, stressed VaR value
#' @param s dbl, stressed CVaR value
#' @param kappa dbl, compound Poisson parameter
#' @param grid_max int
#' @param time_stress dbl
#' @param N nt
#'
#' @return
#' @export
#'
eta2_eq <- function(eta, dist, q, s, kappa, grid_max, time_stress=1, N=8192) {
  # scale variables to avoid hitting numerical max
  terminal_eqn <- function(x) {
    q_scaled <- q / grid_max
    s_scaled <- s / grid_max
    x_scaled <- x / grid_max
    exp(-eta * (x_scaled - q_scaled)) * (x_scaled - s_scaled) *
      as.integer(x_scaled >= q_scaled)
  }
  FST(x=0, t=0,  dist=dist, kappa = kappa, terminal_cond=terminal_eqn,
      grid_max=grid_max, time_stress=time_stress, N=N)
}


#' Computes the optimal eta for CVaR stress
#'
#' @param kappa dbl, compound Poisson parameter
#' @param stress_parms list, c = level of VaR & CVaR, q = stressed value of VaR,
#' s = stressed value of CVaR
#' @param dist RPS_dist object
#' @param time_stress dbl
#' @param N int
#'
#' @return
#' @export
#'
eta_CVaR <- function(kappa,stress_parms,dist,time_stress=1,N=8192) {
  with(stress_parms,{
    grid_max <- 5 * kappa * mean(dist)
    # compute eta 2 using helper function and uniroot
    eta2 <- uniroot(function(eta) eta2_eq(eta, dist=dist, q=q, s=s,
                                          kappa=kappa, grid_max= grid_max,
                                          time_stress=time_stress, N=N),
                    interval = c(-5,5),tol=1e-10, extendInt = "yes")$root
    eta2 <- eta2 / grid_max
    # compute eta 1
    p <- P_prob(x=0,t=0,kappa=kappa,q=q,dist=dist,
                grid_max=grid_max,time_stress=time_stress,N=N)
    expectation <- P_expectation(x=0,t=0,kappa=kappa,q=q,eta2=eta2,dist=dist,
                                 grid_max=grid_max,time_stress=time_stress,N=N)
    eta1 <- log((1-c)/c * p/expectation)

    # return eta1 and eta2
    c(eta1,eta2)
  })
}


#' Computes h(t,x,y) for CVaR constraint
#'
#' @param x dbl
#' @param y dbl
#' @param t dbl
#' @param kappa dbl, compound Poisson parameter
#' @param q dbl, stressed VaR value
#' @param dist RPS_dist objecrt
#' @param eta dbl
#' @param stress_parms list of parameters for stress; must contain:
#' q: stressed VaR value, time_stress: time at which stress occurs
#' @param N int
#'
#' @return
#' @export
#'
h_CVaR <- function(x,y,t,kappa,dist,eta,stress_parms,N=8192){
  with(stress_parms,{
    if (t > time_stress) h <- 1
    else{
      grid_max <- max(x) + max(y) + 30
      num <- exp(-eta[1]) * P_prob(x=x+y,t=t,kappa=kappa,q=q,
                                   dist=dist,grid_max=grid_max,
                                   time_stress=time_stress,N=N) +
        P_expectation(x=x+y,t=t,kappa=kappa,q=q,eta2=eta[2],dist=dist,
                      grid_max=grid_max,time_stress=time_stress,N=N)
      denom <- exp(-eta[1]) * P_prob(x=x,t=t,kappa=kappa,q=q,
                                     dist=dist,grid_max=grid_max,
                                     time_stress=time_stress,N=N) +
        P_expectation(x=x,t=t,kappa=kappa,q=q,eta2=eta[2],dist=dist,
                      grid_max=grid_max,time_stress=time_stress,N=N)
      h <- num / denom
      }
    return(h)
  })
}


