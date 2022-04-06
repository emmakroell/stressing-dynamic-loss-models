# functions to compute correlation along the path and plot it

#' compute correlation CI for Spearman's rho and Kendall's tau
#' based on Bonnett and Wright (2000)
#' @param r correlation point estimate
#' @param n number of observations
#' @param method one of "spearman", "kendall"
#' @param alpha confidence level
#'
#' @return
#' @export
#'
corr_CI <- function(r,n,method,alpha=0.05){
  if (method == "spearman") constant <- (1+r^2/2)/(n-3)
  if (method == "kendall") constant <-0.437/(n-4)
  tanh(atanh(r) + c(-1,1) * sqrt(constant)*stats::qnorm(alpha/2))
}



#' Compute correlation with CI along a path
#'
#' @param X1 vector 1
#' @param X2 vector 2
#' @param time_seq times at which to evaluate vectors
#' @param method one of "pearson", "spearman", "kendall"
#'
#' @return
#' @export
#'
pathwise_corr <- function(X1,X2,time_seq,method){
  estimate <- rep(NA, length(time_seq))
  lwr <- rep(NA, length(time_seq))
  upr <- rep(NA, length(time_seq))

  for (i in 1:length(time_seq)){
    corr_temp <- stats::cor.test(X1[time_seq[i],],X2[time_seq[i],],method=method)
    estimate[i] <- corr_temp$estimate

    if (method == "pearson"){
      lwr[i] <- corr_temp$conf.int[1]
      upr[i] <- corr_temp$conf.int[2]
    } else if (method == "spearman"){
      CI <- corr_CI(estimate[i],length(X2[time_seq[i],]),method=method)
      lwr[i] <- CI[1]
      upr[i] <- CI[2]
    } else if (method == "kendall"){
      CI <- corr_CI(estimate[i],length(X2[time_seq[i],]),method=method)
      lwr[i] <- CI[1]
      upr[i] <- CI[2]
    } else{
      stop("Unknown method")
    }
  }
  dplyr::as_tibble(cbind(estimate,lwr,upr,time_seq=time_seq))
}



#' Create plot comparing pathwise correlation under P and Q
#'
#' @param object RPS_model object
#' @param nsteps number of steps at which to evaluate the correlation
#' @param method one of "pearson", "spearman", "kendall"
#'
#' @return
#' @export
#'
compare_correlation <- function(object,nsteps,method){
  time_seq <- seq(1,length(object$time_vec),length.out=nsteps)
  baseline_sim <- sim_baseline_biv(object)
  X1_P <- baseline_sim$X1
  X2_P <- baseline_sim$X2

  cor1 <- cbind(pathwise_corr(object$paths$X1,object$paths$X2,
                              time_seq,method=method),
                time = object$time_vec[time_seq],
                measure="stressed")
  cor2 <- cbind(pathwise_corr(X1_P,X2_P,time_seq,method=method),
                time = object$time_vec[time_seq],
                measure="original")

  plot <- rbind(cor1,cor2) %>%
    dplyr::mutate(dplyr::across(1:4, as.numeric))%>%
    ggplot2::ggplot(aes(x=time,y=estimate,ymin=lwr,ymax=upr,
                        colour=measure,fill=measure)) +
    ggplot2::geom_line() + ggplot2::geom_ribbon(alpha=0.5) +
    ggplot2::theme_bw()

  return(list(cor = rbind(cor1,cor2), plot = plot))
}


# compare_correlation(gamma_t_ex2,10,method="spearman")
# compare_correlation(gamma_ind_ex2,50,method="spearman")
#
# res1 <- compare_correlation(gamma_ind_ex2,50)
# res2 <- compare_correlation(gamma_t_ex2,50)
#
# headings <- c('t' = "t copula", 'ind' = "Indpedent copula")
#
# rbind(cbind(res1$cor,copula = "ind"),
#       cbind(res2$cor,copula = "t")) %>%
#   mutate(across(1:4, as.numeric))%>%
#   ggplot(aes(x=time,y=estimate,ymin=lwr,ymax=upr)) +
#   geom_line(aes(colour=measure)) + geom_ribbon(aes(fill=measure),alpha=0.5) +
#   facet_wrap(~copula,labeller = ggplot2::as_labeller(headings)) + theme_bw()
#
#
#
# res1 <- compare_correlation(gamma_ind_ex2,50,method="spearman")
# res2 <- compare_correlation(gamma_t_ex2,50,method="spearman")
# res3 <- compare_correlation(gamma_counter_ex,50,method="spearman")
# res4 <- compare_correlation(gamma_comono_ex,50,method="spearman")
#
# headings <- c('t' = "t copula", 'ind' = "Indpedent copula",
#               "counter" = "Countermonotonic","comono" = "Comonotonic")
#
# rbind(cbind(res1$cor,copula = "ind"),
#       cbind(res2$cor,copula = "t"),
#       cbind(res3$cor,copula = "counter"),
#       cbind(res4$cor,copula = "comono")) %>%
#   mutate(across(1:4, as.numeric))%>%
#   ggplot(aes(x=time,y=estimate,ymin=lwr,ymax=upr)) +
#   geom_line(aes(colour=measure)) +
#   geom_ribbon(aes(fill=measure),alpha=0.5) +
#   facet_wrap(~copula,labeller = ggplot2::as_labeller(headings),nrow=1) +
#   theme_bw() + ggtitle(expression("Spearman's" ~ rho))
#
#
# res1 <- compare_correlation(gamma_ind_ex2,20,method="kendall")
# res2 <- compare_correlation(gamma_t_ex2,20,method="kendall")
# res3 <- compare_correlation(gamma_counter_ex,20,method="kendall")
# res4 <- compare_correlation(gamma_comono_ex,20,method="kendall")
#
# headings <- c('t' = "t copula", 'ind' = "Indpedent copula",
#               "counter" = "Countermonotonic","comono" = "Comonotonic")
#
# rbind(cbind(res1$cor,copula = "ind"),
#       cbind(res2$cor,copula = "t"),
#       cbind(res3$cor,copula = "counter"),
#       cbind(res4$cor,copula = "comono")) %>%
#   mutate(across(1:4, as.numeric))%>%
#   ggplot(aes(x=time,y=estimate,ymin=lwr,ymax=upr)) +
#   geom_line(aes(colour=measure)) +
#   geom_ribbon(aes(fill=measure),alpha=0.5) +
#   facet_wrap(~copula,labeller = ggplot2::as_labeller(headings),nrow=1) +
#   theme_bw() + ggtitle(expression("Kendall's" ~ tau))
#

