# compute empirical copula and plot

#' Compute empirical copula
#'
#' @param X1 vector
#' @param X2 vector
#'
#' @return
#' @export
#'
emp_copula <- function(X1,X2){
  ecdf1 <- stats::ecdf(X1)
  U1 <- ecdf1(X1)
  ecdf2 <- stats::ecdf(X2)
  U2 <- ecdf2(X2)

  data.frame(U1 = U1, U2 = U2)
}



# plot_copula <- function(object,time){
#   time_index <- match(time,object$time_vec)
#
#   # stressed
#   X1 <- object$paths$X1[time_index,]
#   X2 <- object$paths$X2[time_index,]
#   stressed <- emp_copula(X1,X2) %>%
#     dplyr::mutate(type = "stressed")
#
#   # baseline
#   baseline <- sim_baseline_biv(object)
#   X1 <- baseline$X1[time_index,]
#   X2 <- baseline$X2[time_index,]
#   baseline <- emp_copula(X1,X2) %>%
#     dplyr::mutate(type = "baseline")
#
#   rbind(baseline,stressed) %>%
#     ggplot2::ggplot(ggplot2::aes(U1,U2,colour=type,alpha=type)) +
#     ggplot2::geom_point(size=1) +
#     ggplot2::theme_bw(base_size=14) +
#     ggplot2::scale_alpha_discrete(range=c(0.5,0.3))
# }


#' Plot empirical copula under P and Q
#'
#' @param object RPS_model object
#' @param time time at which to plot
#'
#' @return
#' @export
#'
#' @examples
plot_copula <- function(object,time){
  time_index <- match(time,object$time_vec)

  # stressed
  X1 <- object$paths$X1[time_index,]
  X2 <- object$paths$X2[time_index,]
  stressed <- emp_copula(X1,X2) %>%
    dplyr::mutate(type = "stressed")

  # baseline
  baseline <- sim_baseline_biv(object)
  X1 <- baseline$X1[time_index,]
  X2 <- baseline$X2[time_index,]
  baseline <- emp_copula(X1,X2) %>%
    dplyr::mutate(type = "baseline")

  rbind(baseline,stressed) %>%
    ggplot2::ggplot(ggplot2::aes(U1,U2,colour=type)) +
    ggplot2::facet_grid(~type)+
    ggplot2::geom_point(size=1,alpha=0.4) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none")
    # ggplot2::scale_alpha_discrete(range=c(0.7,0.4))
}



#
# # independent copula
# cop <- indepCopula(dim=2)
# dist <-  mvdc(copula=cop, margins=c("gamma", "gamma"),
#               paramMargins=list(list(shape=2, scale=1),
#                                 list(shape=2, scale=1)))
# U <- rMvdc(5000, dist)
# df <- emp_copula(U[,1],U[,2])
# p <- ggplot(df, aes(U1, U2)) + geom_point() + theme_classic()
# ggExtra::ggMarginal(p, type = "histogram")
#
#
# cop <- tCopula(0.8,df=3,dim=2)
# dist <-  mvdc(copula=cop, margins=c("gamma", "gamma"),
#               paramMargins=list(list(shape=2, scale=1),
#                                 list(shape=2, scale=1)))
# U <- rMvdc(5000, dist)
# df <- emp_copula(U[,1],U[,2])
# p <- ggplot(df, aes(U1, U2)) + geom_point() + theme_classic()
# ggExtra::ggMarginal(p, type = "histogram")
#
#
# cop <- lowfhCopula(dim=2)
# dist <-  mvdc(copula=cop, margins=c("gamma", "gamma"),
#               paramMargins=list(list(shape=2, scale=1),
#                                 list(shape=2, scale=1)))
# U <- rMvdc(5000, dist)
# df <- emp_copula(U[,1],U[,2])
# p <- ggplot(df, aes(U1, U2)) + geom_point() + theme_classic()
# ggExtra::ggMarginal(p, type = "histogram")


# plot_copula(gamma_ind_ex1,1)
#
# plot_copula(gamma_t_2_ex3,1)
#
# ggpubr::ggarrange(plot_copula(gamma_ind_ex1,1),plot_copula(gamma_t_2_ex3,1),
#                   common.legend = TRUE)
#
# baseline <- sim_baseline_biv(gamma_ind_ex1)
#
# X1 <- baseline$X1[501,]
# X2 <- baseline$X2[501,]
#
# ecdf1 <- stats::ecdf(X1)
# U1 <- ecdf1(X1)
#
# ecdf2 <- stats::ecdf(X2)
# U2 <- ecdf2(X2)
#
# data.frame(U1 = U1, U2 = U2) %>%
#   ggplot2::ggplot(ggplot2::aes(U1, U2)) +
#   ggplot2::geom_point() +
#   ggplot2::theme_bw()
#
# plot_copula(gamma_ind_ex1,1)
# plot_copula(gamma_ind_ex2,1)
# plot_copula(gamma_ind_ex3,1)
# plot_copula(gamma_t_ex1,1)
# plot_copula(gamma_t_ex2,1)
# plot_copula(gamma_t_ex3,1)
#
#
#
#
#
#
