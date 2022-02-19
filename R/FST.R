
#' Fourier Space Time-stepping (FST) algorithm
#'
#' adapted from the Fourier Space Time-stepping (FST) library,
#' Copyright (C) 2007  Vladimir Surkov
#' https://sourceforge.net/projects/fst-framework/
#' license is GPL-3
#'
#' @param x vector, where to evaluate the FST
#' @param t double, time at which to evaluate the FST
#' @param dist RPS_dist object, contains jump size distribution
#' @param kappa double, compound Poisson parameter
#' @param terminal_cond double, end value
#' @param grid_max double, max grid value for FST
#' @param endtime double, end time
#' @param N int, number of points in grid
#' @param h Boolean, function returns fraction of FST at x+y/FST at x if TRUE
#' @param y vector, optional argument for numerator if h = TRUE
#'
#' @return FST approximation of terminal_cond at x at time t
#' @export
#'
FST <- function(x, t, dist, kappa, terminal_cond, grid_max, endtime=1, N=5e3,
                h = FALSE, y = NULL){
  with(dist, {
    # Real space
    x_min <- 0
    x_max <- grid_max

    dx <- (x_max-x_min)/(N-1)
    x.seq <- seq(x_min,x_max,by=dx)

    # Fourier space
    w_max <- pi/dx
    dw <- 2*w_max/N
    w <- c(seq(0,w_max,by=dw),seq(-w_max+dw,-dw,by=dw))

    # terminal condition
    f <- terminal_cond(x.seq)

    # FST method
    char_exp <- kappa * (char_fun(w,parms) - 1)
    char_exp_factor <-  exp(char_exp * (endtime-t))

    # compute omega
    omega <- Re(fft(fft(f)*char_exp_factor, inverse = TRUE)/length(f))

    # interpolate
    if (h == TRUE) {
      num <- spline(x=x.seq,y=omega,xout=x+y)$y
      denom <- spline(x=x.seq,y=omega,xout=x)$y
      res <- num/denom
    } else {
      res <- spline(x=x.seq,y=omega,xout=x)$y
    }
    return(res)
  })
}


