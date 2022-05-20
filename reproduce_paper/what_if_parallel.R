library(doParallel)
library(tidyverse)
library(copula)
devtools::load_all()
source("R/class.R")
source("R/FST.R")
source("R/FST_implementations.R")
source("R/2d_sim_model.R")
source("R/run_parallel.R")

# set parms
kappa <- 5
beta <- 0.9
perc_target <- 5
dt <- 2e-3
tol <- 5e-4

# jump distribution
gamma_exp_t <- new_RPS_dist_biv(copula=tCopula(0.8,df=3,dim=2),
                                margins=c("gamma", "gamma"),
                                sim_fun = function(x,parms) rgamma(x,shape=parms$alpha1,
                                                                   rate=parms$beta1),
                                dens_fun = function(x,parms) dgamma(x,shape=parms$alpha1,
                                                                    rate=parms$beta1),
                                char_fun = function(x,parms) (1 - 1i*x/parms$beta1)^(-parms$alpha1),
                                mean_fun = function(parms) parms$alpha1/parms$beta1,
                                parms = list(alpha1 = 2, beta1=1, alpha2=1, beta2=2))

# compute baseline VaR sum
baseline <- sim_baseline_biv(jump_dist=gamma_exp_t,kappa=5,
                             Npaths=1e4,time_vec=seq(0,1,by=dt))
sum_VaR_P <- quantile(baseline$X1[length(seq(0,1,by=dt)),] + baseline$X2[length(seq(0,1,by=dt)),],
                      beta,type=1)


# ------------------------------------------------------------------------------

# find stress - alpha = 0.4
res4 <- find_stress_parallel(floor=1.25, ceiling=1.55,
                             perc_target=perc_target,
                             sum_VaR_P=as.numeric(sum_VaR_P),
                             alpha=0.4, beta=beta,
                             jump_dist=gamma_exp_t, kappa=kappa,
                             dt=dt, tol=tol, max_iter=4)
saveRDS(res4,"saved_sims/stress_find/stress_alpha_04.RDS")

# find stress - alpha = 0.5
res5 <- find_stress_parallel(floor=1.2, ceiling=1.5,
                             perc_target=perc_target,
                             sum_VaR_P=as.numeric(sum_VaR_P),
                             alpha=0.5, beta=beta,
                             jump_dist=gamma_exp_t, kappa=kappa,
                             dt=dt, tol=tol, max_iter=4)
saveRDS(res5,"saved_sims/stress_find/stress_alpha_05.RDS")

# find stress - alpha = 0.6
res6 <- find_stress_parallel(floor=1.15, ceiling=1.45,
                             perc_target=perc_target,
                             sum_VaR_P=as.numeric(sum_VaR_P),
                             alpha=0.6, beta=beta,
                             jump_dist=gamma_exp_t, kappa=kappa,
                             dt=dt, tol=tol, max_iter=4)
saveRDS(res6,"saved_sims/stress_find/stress_alpha_06.RDS")

# find stress - alpha = 0.7
res7 <- find_stress_parallel(floor=1.1, ceiling=1.4,
                             perc_target=perc_target,
                             sum_VaR_P=as.numeric(sum_VaR_P),
                             alpha=0.7, beta=beta,
                             jump_dist=gamma_exp_t, kappa=kappa,
                             dt=dt, tol=tol, max_iter=4)
saveRDS(res7,"saved_sims/stress_find/stress_alpha_07.RDS")

# find stress - alpha = 0.8
res8 <- find_stress_parallel(floor=1.05, ceiling=1.35,
                             perc_target=perc_target,
                             sum_VaR_P=as.numeric(sum_VaR_P),
                             alpha=0.8, beta=beta,
                             jump_dist=gamma_exp_t, kappa=kappa,
                             dt=dt, tol=tol, max_iter=4)
saveRDS(res8,"saved_sims/stress_find/stress_alpha_08.RDS")

# find stress - alpha = 0.9
res9 <- find_stress_parallel(floor=1.05, ceiling=1.35,
                             perc_target=perc_target,
                             sum_VaR_P=as.numeric(sum_VaR_P),
                             alpha=0.9, beta=beta,
                             jump_dist=gamma_exp_t, kappa=kappa,
                             dt=dt, tol=tol, max_iter=4)
saveRDS(res9,"saved_sims/stress_find/stress_alpha_09.RDS")


# find stress - alpha = 0.3
res3 <- find_stress_parallel(floor=1.4, ceiling=1.7,
                             perc_target=perc_target,
                             sum_VaR_P=as.numeric(sum_VaR_P),
                             alpha=0.3, beta=beta,
                             jump_dist=gamma_exp_t, kappa=kappa,
                             dt=dt, tol=tol, max_iter=4)
saveRDS(res3,"saved_sims/stress_find/stress_alpha_03.RDS")

# find stress - alpha = 0.2
res2 <- find_stress_parallel(floor=1.5, ceiling=1.8,
                             perc_target=perc_target,
                             sum_VaR_P=as.numeric(sum_VaR_P),
                             alpha=0.2, beta=beta,
                             jump_dist=gamma_exp_t, kappa=kappa,
                             dt=dt, tol=tol, max_iter=4)
saveRDS(res2,"saved_sims/stress_find/stress_alpha_02.RDS")

# find stress - alpha = 0.1
res1 <- find_stress_parallel(floor=1.55, ceiling=1.85,
                             perc_target=perc_target,
                             sum_VaR_P=as.numeric(sum_VaR_P),
                             alpha=0.1, beta=beta,
                             jump_dist=gamma_exp_t, kappa=kappa,
                             dt=dt, tol=tol, max_iter=4)
saveRDS(res1,"saved_sims/stress_find/stress_alpha_01.RDS")


# ------------------------------------------------------------------------------
# make table with results
res6$table %>%
  mutate(diff = abs(perc_incr - perc_target)) %>%
  filter(diff == min(diff)) %>% summarise_all(,.funs=mean)

res4$table %>% mutate(diff = abs(perc_incr - perc_target)) %>%
  filter(diff == min(diff))
res4$table %>% mutate(diff = abs(perc_incr - perc_target)) %>%
  filter(diff == min(diff)) %>% select(stress) %>%  pull()
# only good to 2 decimal places
res4_min <- res4$table %>%
  mutate(diff = abs(perc_incr - perc_target)) %>%
  filter(diff == min(diff))

res3_min <- res3$table %>%
  mutate(diff = abs(perc_incr - perc_target)) %>%
  filter(diff == min(diff))

res5_min <- res5$table %>%
  mutate(diff = abs(perc_incr - perc_target)) %>%
  filter(diff == min(diff))

res6_min <- res6$table %>%
  mutate(diff = abs(perc_incr - perc_target)) %>%
  filter(diff == min(diff))

res7_min <- res7$table %>%
  mutate(diff = abs(perc_incr - perc_target)) %>%
  filter(diff == min(diff))

res8_min <- res8$table %>%
  mutate(diff = abs(perc_incr - perc_target)) %>%
  filter(diff == min(diff))

res9_min <- res9$table %>%
  mutate(diff = abs(perc_incr - perc_target)) %>%
  filter(diff == min(diff))

results_table <- rbind(res3_min[1,],res4_min[1,],res5_min[1,],res6_min[1,],
                       res7_min[1,],res8_min[1,],res9_min[1,])

results_table %>% mutate(sum_VaR_Q = sum_VaR) %>%
  select(alpha,stress,sum_VaR_P,sum_VaR_Q,perc_incr,diff)

results_table %>% mutate(sum_VaR_Q = sum_VaR) %>%
  select(alpha,stress,sum_VaR_P,sum_VaR_Q,perc_incr,diff) %>%
  mutate(across(2:6, function(x) round(x,2))) %>%
  kableExtra::kbl() %>% kableExtra::kable_styling()
  # kableExtra::kbl("latex")
