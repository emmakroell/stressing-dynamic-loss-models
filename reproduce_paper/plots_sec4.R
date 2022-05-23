# Plots for section 4 of paper
devtools::load_all()

gamma_2_1 <- new_RPS_dist(dist_fun = function(x,parms) pgamma(x,shape=parms$alpha,
                                                              rate=parms$beta),
                          dens_fun = function(x,parms) dgamma(x,shape=parms$alpha,
                                                              rate=parms$beta),
                          sim_fun = function(x,parms) rgamma(x,shape=parms$alpha,
                                                             rate=parms$beta),
                          char_fun = function(x,parms) (1 - 1i*x/parms$beta)^(-parms$alpha),
                          min = 0,
                          max = qgamma(1-1e-10,shape=2,rate=1),
                          mean_fun = function(parms) parms$alpha/parms$beta,
                          parms = list(alpha = 2, beta=1))

# ex_gamma_VaR <- stressed_sim(kappa = 5, jump_dist = gamma_2_1,
#                              stress_type = "VaR",
#                              stress_parms = list(c=0.9,VaR_stress=1.1),
#                              Npaths=1e4, end_time=1, dt=1e-3)

ex_gamma_VaR <- readRDS("saved_sims/1d_VaR_10.RDS")

ex_gamma_VaR$stress_parms$q
quantile(ex_gamma_VaR$paths[length(ex_gamma_VaR$time_vec),],0.9,type=1)
baseline_VaR <- sim_baseline(ex_gamma_VaR)
quantile(baseline_VaR[length(ex_gamma_VaR$time_vec),],0.9,type=1)
FNN::KL.divergence(ex_gamma_VaR$paths[length(ex_gamma_VaR$time_vec),],
                   baseline_VaR[length(ex_gamma_VaR$time_vec),],100)

# saveRDS(ex_gamma_VaR, file = "1d_VaR_10.RDS")

# find useful paths
get_paths <- function(object,cutpoints){
  paths_to_plot <- rep(NA,length(cutpoints)-1)
  withr::with_seed(4529,{
    for (i in 2:length(cutpoints)){
      paths_to_plot[i-1] <- sample(which(object$paths[length(object$time_vec),] < cutpoints[i] &
                                         object$paths[length(object$time_vec),] > cutpoints[i-1]),1)
    }
  })
  return(paths_to_plot)
}

my_paths <- get_paths(ex_gamma_VaR,c(0,4,10,18,19,25))
#
# plot_paths(ex_gamma_VaR, Npaths=5,
#            quantiles=list(lower=0.1,upper=0.9),
#            indices=my_paths)
# pdf 5 by 8
plot_paths2(ex_gamma_VaR, Npaths=5,
            quantiles=list(lower=0.1,upper=0.9),
            indices=my_paths)

plot_G_Q(x_vec=c(16,18,20),t_vec=c(0.5,1),kappa=5,
         jump_dist=gamma_2_1,stress_type="VaR",
         stress_parms = list(c=0.9,VaR_stress=1.1),xmax=8)
# pdf 5 by 8


# ex_gamma_CVaR <- stressed_sim(kappa = 5, jump_dist = gamma_2_1,
#                               stress_type = "CVaR",
#                               stress_parms = list(c=0.9,
#                                                   VaR_stress=1.1,
#                                                   CVaR_stress=1.08,
#                                                   time_stress = 1),
#                               Npaths=1e4, end_time=1, dt=1e-3)

# saveRDS(ex_gamma_CVaR, file = "1d_CVaR_10.RDS")

ex_gamma_CVaR <- readRDS("saved_sims/1d_CVaR_10.RDS")

plot_paths2(ex_gamma_CVaR, Npaths=5,
           quantiles=list(lower=0.1,upper=0.9),
           indices=my_paths)
plot_G_Q(x_vec=c(16,18,20),t_vec=c(0.5,1),kappa=5,
         jump_dist=gamma_2_1,stress_type="CVaR",
         stress_parms = list(c=0.9,VaR_stress=1.1,CVaR_stress=1.08),xmax=8)

ex_gamma_CVaR$stress_parms$q
quantile(ex_gamma_CVaR$paths[length(ex_gamma_CVaR$time_vec),],0.9,type=1)
baseline_CVaR <- sim_baseline(ex_gamma_CVaR)
quantile(baseline_CVaR[length(ex_gamma_CVaR$time_vec),],0.9,type=1)
FNN::KL.divergence(ex_gamma_CVaR$paths[length(ex_gamma_CVaR$time_vec),],
                   baseline_CVaR[length(ex_gamma_CVaR$time_vec),],100)

