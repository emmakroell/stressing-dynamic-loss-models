library(tidyverse)
library(copula)
devtools::load_all()

gamma_ind <- new_RPS_dist_biv(copula=indepCopula(dim=2),
                              margins=c("gamma", "gamma"),
                              sim_fun = function(x,parms) rgamma(x,shape=parms$alpha1,
                                                                 rate=parms$beta1),
                              dens_fun = function(x,parms) dgamma(x,shape=parms$alpha1,
                                                                  rate=parms$beta1),
                              char_fun = function(x,parms) (1 - 1i*x/parms$beta1)^(-parms$alpha1),
                              mean_fun = function(parms) parms$alpha1/parms$beta1,
                              parms = list(alpha1 = 2, beta1=1, alpha2=1, beta2=2))


gamma_exp_t <- new_RPS_dist_biv(copula=tCopula(0.8,df=3,dim=2),
                                margins=c("gamma", "gamma"),
                                sim_fun = function(x,parms) rgamma(x,shape=parms$alpha1,
                                                                   rate=parms$beta1),
                                dens_fun = function(x,parms) dgamma(x,shape=parms$alpha1,
                                                                    rate=parms$beta1),
                                char_fun = function(x,parms) (1 - 1i*x/parms$beta1)^(-parms$alpha1),
                                mean_fun = function(parms) parms$alpha1/parms$beta1,
                                parms = list(alpha1 = 2, beta1=1, alpha2=1, beta2=2))


# using plotting functions
plot1 <- plot_G_marginals(
  x_vec = c(16,18,20),
  t_vec = 1,
  kappa = 5,
  jump_dist = gamma_ind,
  stress_type = "VaR",
  stress_parms = list(c=0.9,VaR_stress=1.15,time_stress=1),
  N_out = 10000,
  xmax = c(9,4)
)$plot_y2
plot2 <- plot_G_marginals(
  x_vec = c(16,18,20),
  t_vec = 1,
  kappa = 5,
  jump_dist = gamma_exp_t,
  stress_type = "VaR",
  stress_parms = list(c=0.9,VaR_stress=1.15,time_stress=1),
  N_out = 10000,
  xmax = c(9,4)
)$plot_y2
ggpubr::ggarrange(plot1,plot2)


# plot for paper
ind_cop <- plot_G_marginals(
  x_vec = c(16,18,20),
  t_vec = 1,
  kappa = 5,
  jump_dist = gamma_ind,
  stress_type = "VaR",
  stress_parms = list(c=0.9,VaR_stress=1.15,time_stress=1),
  N_out = 10000,
  xmax = c(9,3))$plot_data

t_cop <- plot_G_marginals(
  x_vec = c(16,18,20),
  t_vec = 1,
  kappa = 5,
  jump_dist = gamma_exp_t,
  stress_type = "VaR",
  stress_parms = list(c=0.9,VaR_stress=1.15,time_stress=1),
  N_out = 10000,
  xmax = c(9,3))$plot_data

x_P <- seq(0.09,3,0.05)
y_P <- dgamma(x, shape=1, rate=2)

rbind(cbind(ind_cop,copula="indep"),
      cbind(t_cop,copula="t")) %>%
  dplyr::filter(dim == 2, time == 1) %>%
  ggplot2::ggplot() +
  ggplot2::geom_histogram(ggplot2::aes(x=value,y=..density..,colour=x,fill=x),
                          bins=25,alpha=0.6) +
  ggplot2::geom_line(data=dplyr::as_tibble(cbind(x_P,y_P)),
                     ggplot2::aes(x=x_P,y=y_P)) +
  ggplot2::facet_grid(copula ~ x,
                      labeller = as_labeller(c("16"="x[1] == 16",
                                               "18"="x[1] == 18",
                                               "20"="x[1] == 20",
                                               "indep" = "Independence ~ copula",
                                               "t" = "t-copula"),
                                             ggplot2::label_parsed)) +
  ggplot2::theme_bw(base_size = 16) +
  ggplot2::theme(legend.position = "none") +
  ggplot2::xlab('') + ggplot2::ylab('') + ggplot2::xlim(c(0,3))

# save as pdf 5 by 8



