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
#' @param endtime dbl
#' @param N int
#'
#' @return
#' @export
#'
P_prob <- function(x,t,kappa,q,dist,grid_max=50,endtime=1,N=8192){
  indic_func <- function(x) as.integer(x < q)
  FST(x=x, t=t, dist=dist, kappa = kappa,terminal_cond=indic_func,
      grid_max=grid_max, endtime=endtime, N=N)
}


#' computes VaR under the P measure using FST and stats::uniroot
#'
#' @param kappa dbl, compound Poisson parameter
#' @param c dbl, level for VaR between 0 and 1
#' @param dist RPS_dist object
#'
#' @return
#' @export
#'
compute_VaR <- function(kappa,c,dist){
  grid_max <- 5 * kappa * mean(dist)
  stats::uniroot(function(q) P_prob(x=0,t=0,kappa=kappa,q=q,
                                    dist=dist,grid_max=grid_max)-c,
                 interval=c(0, 5*kappa*mean(dist)))$root
}


#' Computes the optimal eta for VaR stress using FST
#'
#' @param kappa double, compound Poisson parameter
#' @param stress_parms list, c = level of VaR, q = stressed value of VaR
#' @param dist RPS_dist object, contains jump size distribution
#'
#' @return double, optimal eta for VaR stress
#' @export
#'
eta_VaR <- function(kappa,stress_parms,dist) {
  with(stress_parms,{
    grid_max <- 5 * kappa * mean(dist)
    p <- P_prob(x=0, t=0, kappa = kappa, q=q, dist=dist, grid_max=grid_max)
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
#' @param c double, VaR level to stress
#' @param kappa double, compound Poisson parameter
#' @param endtime double, end time
#' @param N int, number of points for FST method
#'
#' @return
#' @export
#'
h_VaR <- function(x, y, t, dist, eta, q, c, kappa, endtime=1, N=8192) {
  grid_max <- max(x) + max(y) + 30
  num <- 1 + (exp(-eta) - 1) * P_prob(x+y,t,kappa,q,dist,grid_max,endtime,N)
  denom <- 1 + (exp(-eta) - 1) * P_prob(x,t,kappa,q,dist,grid_max,endtime,N)
  num / denom
}



