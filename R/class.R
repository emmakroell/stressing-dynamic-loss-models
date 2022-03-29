#' Defines RPS_model class
#'
#' @param jump_dist RPS_dist object, contains jump size distribution
#' @param kappa double, compound Poisson parameter
#' @param stress_type string, type of stress, one of "VaR", "VaR & CVaR"
#' @param stress_parms list, parameters for stress
#' @param paths matrix, paths of stressed model
#' @param time_vec vector, times at which stressed model is evaluated
#' @param kappa_Q matrix, paths of stressed processed intensity
#'
#' @return
#' @export
#'
new_RPS_model <- function(jump_dist, kappa, stress_type,
                          stress_parms, paths, time_vec, kappa_Q){

  stress_type <- match.arg(stress_type, c("VaR", "CVaR"))

  model <- list(jump_dist = jump_dist,
                kappa = kappa,
                stress_type = stress_type, # type of stress: VaR or VaR & CVaR
                stress_parms = stress_parms, # list of stress parameters corresponding to stress type
                paths = paths,
                time_vec = time_vec,
                kappa_Q = kappa_Q
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
new_RPS_dist <- function(copula, margins, dens_fun, char_fun, mean_fun, parms) {
  # copula <- copula(param = NULL, dim=2) #copula(0.5,dim=2) #for normal FIX
  biv_dist = mvdc(copula=copula,
                  margins=margins,
                  paramMargins=list(list(shape=parms$alpha1,
                                         scale=parms$beta1),
                                    list(shape=parms$alpha2,
                                         scale=parms$beta2)))

  model <- list(biv_dist = biv_dist,
                dens_fun = dens_fun,
                char_fun = char_fun,
                mean_fun = mean_fun,
                parms = parms
  )
  ## Name of the class
  attr(model, "class") <- "RPS_dist"
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
is.RPS_dist <- function(object) inherits(object, "RPS_dist")

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
