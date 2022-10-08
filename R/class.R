#' Defines RPS_model class
#'
#' @param jump_dist RPS_dist object, contains jump size distribution
#' @param kappa double, compound Poisson parameter
#' @param stress_type string, type of stress, one of "VaR", "VaR & CVaR"
#' @param stress_parms list, parameters for stress
#' VaR: q: stressed VaR value, c: VaR level, VaR_stress: multiplier for P-VaR
#' CVAR: q: stressed VaR value, c: VaR level, s: stressed CVaR level,
#' VaR_stress: multiplier for P-VaR, CVaR_stress: multiplier for P-CVaR
#' time_stress: time at which stress occurs; if NULL, assumed to be end_time
#' @param paths matrix, paths of stressed model
#' @param time_vec vector, times at which stressed model is evaluated
#' @param intensity matrix, paths of stressed processed intensity
#' @param model_type string, one of "univariate", "bivariate_copula", or
#'                   "bivariate_mixture"
#'
#' @return
#' @export
#'
RPS_model <- function(model_type, jump_dist, kappa, stress_type, eta,
                          stress_parms, paths, time_vec, intensity){

  stress_type <- match.arg(stress_type, c("VaR", "CVaR"))

  model <- list(model_type = model_type,
                jump_dist = jump_dist,
                kappa = kappa,
                stress_type = stress_type, # type of stress: VaR or VaR & CVaR
                stress_parms = stress_parms, # list of stress parameters corresponding to stress type
                paths = paths,
                eta = eta,
                time_vec = time_vec,
                intensity = intensity
  )
  ## Name of the class
  attr(model, "class") <- "RPS_model"
  return(model)
}

#' Method for RPS_model for is()
#'
#' @param object RPS_model
#'
#' @return Boolean
#' @export
#'
is.RPS_model <- function(object) inherits(object, "RPS_model")


# UNIVARIATE
#' Defines RPS_dist class, class giving a distribution
#'
#' @param dist_fun function of x and parms, returns distribution function
#' @param dens_fun function of x and parms, returns density function
#' @param sim_fun function of x and parms, returns random draw generation function
#' @param char_fun function of x and parms, returns characteristic function
#' @param mean_fun function of parms, returns mean
#' @param parms list, parameters for distribution
#' @param min dbl, where to start integration of random variable
#' @param max dbl, where to end integration of random variable
#'
#' @return
#' @export
#'
RPS_dist_univ <- function(dist_fun, dens_fun, sim_fun, char_fun, mean_fun,
                              min, max, parms) {

  model <- list(dist_fun = dist_fun,
                dens_fun = dens_fun,
                sim_fun = sim_fun,
                char_fun = char_fun,
                mean_fun = mean_fun,
                min = min,  # approx of ess inf
                max = max,  # approx of ess sup
                parms = parms
  )
  ## Name of the class
  attr(model, "class") <- "RPS_dist_univ"
  return(model)
}

#' Method for RPS_dist for is()
#'
#' @param object RPS_dist
#'
#' @return Boolean
#' @export
#'
is.RPS_dist_univ <- function(object) inherits(object, "RPS_dist_univ")

#' Method for RPS_dist for mean()
#'
#' @param object RPS_dist
#'
#' @return mean of distribution
#' @export
#'
#' @example
#' new_RPS_dist(dist_fun = function(x,parms) pgamma(x,shape=parms$alpha,rate=parms$beta),
#'              dens_fun = function(x,parms) dgamma(x,shape=parms$alpha,rate=parms$beta),
#'              sim_fun = function(x,parms) rgamma(x,shape=parms$alpha,rate=parms$beta),
#'              char_fun = function(x,parms) (1 - 1i*x/parms$beta)^(-parms$alpha),
#'              mean_fun = function(parms) parms$alpha/parms$beta,
#'              parms = list(alpha = 2, beta=1))
#'
mean.RPS_dist <- function(object) {
  with(object, mean_fun(parms))
}


# BIVARIATE

#' Defines RPS_dist class, class giving a distribution
#'
#' @param dist_fun function of x and parms, returns distribution function
#' @param dens_fun function of x and parms, returns density function
#' @param sim_fun function of x and parms, returns random draw generation function
#' @param char_fun function of x and parms, returns characteristic function
#' @param mean_fun function of parms, returns mean
#' @param parms list, parameters for distribution
#'
#' @return
#' @export
#'
RPS_dist_biv <- function(copula, margins, sim_fun, dens_fun,
                         char_fun, mean_fun, parms) {

  # build distribution using copula package
  biv_dist = copula::mvdc(copula=copula,
                  margins=margins,
                  paramMargins=list(list(shape = parms$alpha1,
                                         rate  = parms$beta1),
                                    list(shape = parms$alpha2,
                                         rate  = parms$beta2)))

  model <- list(biv_dist = biv_dist,
                sim_fun = sim_fun,
                dens_fun = dens_fun,
                char_fun = char_fun,
                mean_fun = mean_fun,
                parms = parms
  )
  ## Name of the class
  attr(model, "class") <- "RPS_dist_biv"
  return(model)
}

# new_RPS_dist <- function(dist_fun, dens_fun, sim_fun, char_fun, mean_fun, parms) {
#
#   model <- list(dist_fun = dist_fun,
#                 dens_fun = dens_fun,
#                 sim_fun = sim_fun,
#                 char_fun = char_fun,
#                 mean_fun = mean_fun,
#                 parms = parms
#   )
#   ## Name of the class
#   attr(model, "class") <- "RPS_dist"
#   return(model)
# }

#' Method for RPS_dist for is()
#'
#' @param object RPS_dist
#'
#' @return Boolean
#' @export
#'
is.RPS_dist_biv<- function(object) inherits(object, "RPS_dist_biv")

#' Method for RPS_dist for mean()
#'
#' @param object RPS_dist
#'
#' @return mean of distribution
#' @export
#'
#' @example
#' new_RPS_dist(dist_fun = function(x,parms) pgamma(x,shape=parms$alpha,rate=parms$beta),
#'              dens_fun = function(x,parms) dgamma(x,shape=parms$alpha,rate=parms$beta),
#'              sim_fun = function(x,parms) rgamma(x,shape=parms$alpha,rate=parms$beta),
#'              char_fun = function(x,parms) (1 - 1i*x/parms$beta)^(-parms$alpha),
#'              mean_fun = function(parms) parms$alpha/parms$beta,
#'              parms = list(alpha = 2, beta=1))
#'
mean.RPS_dist_biv<- function(object) {
  with(object, mean_fun(parms))
}


## MIXTURE
# Mixture model code
RPS_dist_mix <- function(sim_fun_a, sim_fun_b, dens_fun_a, dens_fun_b,
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

