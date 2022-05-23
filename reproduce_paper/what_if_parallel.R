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
res8 <- find_stress_parallel(floor=1.1, ceiling=1.4,
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

# ------------------------------------------------------------------------------

# find stress - alpha = 0.8
res8b <- find_stress_parallel(floor=1.12, ceiling=1.22,
                              perc_target=perc_target,
                              sum_VaR_P=as.numeric(sum_VaR_P),
                              alpha=0.8, beta=beta,
                              jump_dist=gamma_exp_t, kappa=kappa,
                              dt=dt, tol=tol, max_iter=4)
saveRDS(res8b,"saved_sims/stress_find/dt1e3/stress_alpha_08.RDS")

# find stress - alpha = 0.9
res9b <- find_stress_parallel(floor=1.1, ceiling=1.2,
                              perc_target=perc_target,
                              sum_VaR_P=as.numeric(sum_VaR_P),
                              alpha=0.9, beta=beta,
                              jump_dist=gamma_exp_t, kappa=kappa,
                              dt=dt, tol=tol, max_iter=4)
saveRDS(res9b,"saved_sims/stress_find/dt1e3/stress_alpha_09.RDS")


# ------------------------------------------------------------------------------
# make table with results

get_estimate_row <- function(res,perc_target) {
  final_res <- res$table[(6*(res$niter-1)+1):(6*res$niter),]
  floor_ind <- tail(which(final_res$perc_incr < perc_target),1)
  ceiling_ind <- head(which(final_res$perc_incr >= perc_target),1)

  final_res[c(floor_ind,ceiling_ind),] %>% summarise_all(.funs=mean)
}

res_table <- rbind(get_estimate_row(res3,perc_target),
                   get_estimate_row(res4,perc_target),
                   get_estimate_row(res5,perc_target),
                   get_estimate_row(res6,perc_target),
                   get_estimate_row(res7,perc_target),
                   get_estimate_row(res8,perc_target))

res_table$stress

res_table %>% select(alpha,stress) %>%
  mutate(stress = round(stress,3)) %>%
  kableExtra::kbl("latex")


res_table %>% select(alpha,stress) %>%
  mutate(stress = round((stress - 1) * 100,2) )  %>%
  kableExtra::kbl("latex")


