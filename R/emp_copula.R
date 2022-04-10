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


#' Plot empirical copula contour plot under P and Q
#'
#' @param object RPS_model object
#' @param time time at which to plot
#'
#' @return
#' @export
#'
#' @examples
plot_copula_contour <- function(object,time){

  time_index <- match(time,object$time_vec)
  u <- seq(0, 1, length.out = 25)
  grid <- as.matrix(expand.grid("x" = u, "y" = u))


  # stressed
  X1 <- object$paths$X1[time_index,]
  X2 <- object$paths$X2[time_index,]
  draws <- cbind(X1,X2)
  ec <-  copula::C.n(grid, X = draws)
  # contourplot2(cbind(as.matrix(grid),ec))
  res_Q <- cbind(as.matrix(grid),ec)
  colnames(res_Q) <- c("x", "y", "z")
  res_Q <- dplyr::as_tibble(cbind(res_Q,measure="stressed")) %>%
    dplyr::mutate(dplyr::across(1:3,as.numeric))
  plot1 <- ggplot2::ggplot(res_Q, ggplot2::aes(x,y,z=z)) +
    ggplot2::geom_contour_filled() +
    ggplot2::theme_minimal()


  # baseline
  baseline <- sim_baseline_biv(object)
  X1 <- baseline$X1[time_index,]
  X2 <- baseline$X2[time_index,]
  draws <- cbind(X1,X2)
  ec <- copula::C.n(grid, X = draws)
  # contourplot2(cbind(as.matrix(grid),ec))
  res_P <- cbind(as.matrix(grid),ec)
  colnames(res_P) <- c("x", "y", "z")
  res_P <- dplyr::as_tibble(cbind(res_P,measure="original")) %>%
    dplyr::mutate(dplyr::across(1:3,as.numeric))
  plot2 <- ggplot2::ggplot(res_P, ggplot2::aes(x,y,z=z)) +
    ggplot2::geom_contour_filled() +
    ggplot2::theme_minimal()

  xy <- NULL
  for (constant in seq(0.1,0.9,by=0.1)){
    xy <- rbind(xy, as_tibble(list(x = seq(0,1,0.01),
                                   y = sapply(seq(0,1,0.01),
                                              function(x) min(constant/x,1)),
                                   constant = constant)))
  }
  xy <- mutate(xy,across(3,as.numeric),
               y = ifelse(y==1,NA,y))

  # combined plot
  data <- rbind(res_P,res_Q) %>% mutate(measure = fct_rev(measure))
  plot <- data %>%
    ggplot2::ggplot() +
    ggplot2::geom_contour(ggplot2::aes(x,y,z=z,colour=measure), lwd=1) +
    ggplot2::theme_minimal() +
    ggplot2::geom_line(data = xy, aes(x,y,group=constant), lty = 2)

  return(list(data = data,
              plot = plot))
}

#' Plot empirical copula density contour plot under P and Q
#'
#' @param object RPS_model object
#' @param time time at which to plot
#'
#' @return
#' @export
#'
#' @examples
plot_copula_dens_contour <- function(object,time){

  time_index <- match(time,object$time_vec)

  # stressed
  X1 <- object$paths$X1[time_index,]
  X2 <- object$paths$X2[time_index,]
  u1 <- copula::pobs(cbind(X1,X2))
  kde.fit1 <- kdecopula::kdecop(u1)
  plot(kde.fit1)
  graphics::contour(kde.fit1,margins="unif",levels = c(0.1,0.5,1.5,2,3),
                    col="#F8766D",lwd=2,cex.lab=1.5,cex.axis=1.5,
                    labcex=1.1,vfont=c("sans serif", "bold"))

  # baseline
  baseline <- sim_baseline_biv(object)
  X1 <- baseline$X1[time_index,]
  X2 <- baseline$X2[time_index,]
  u2 <- copula::pobs(cbind(X1,X2))
  kde.fit2 <- kdecopula::kdecop(u2)
  plot(kde.fit2)
  graphics::contour(kde.fit2,margins="unif",levels = c(0.1,0.5,1.5,2,3),
                    col="#00BFC4",add=TRUE,lwd=2,
                    labcex=1.1,vfont=c("sans serif", "bold"))
  # add legend
  graphics::legend("bottomright", legend = c("Q", "P"),
         col = c("#F8766D","#00BFC4"),lty=1:1,lwd=c(2,2),cex=1.2)

}


# plot_copula_dens_contour_v2 <- function(object,time){
#
#   time_index <- match(time,object$time_vec)
#
#   # stressed
#   X1 <- object$paths$X1[time_index,]
#   X2 <- object$paths$X2[time_index,]
#   u1 <- copula::pobs(cbind(X1,X2))
#   kde.fit1 <- kdecop(u1)
#
#   contour_dat_Q <- reshape2::melt(kde.fit1$estimate) %>%
#     mutate(u1 = kde.fit1$grid[Var1],
#            u2 = kde.fit1$grid[Var2]) %>%
#     select(u1,u2,value)
#
#   plot1 <- ggplot(contour_dat_Q,aes(x=u1,y=u2,z=value)) +
#     geom_contour(aes(colour = after_stat(level)),
#                  breaks = c(0,0.5,1,1.5,2,3,4,5,10,15))+
#     scale_colour_gradient(low = "red", high = "mediumorchid1") + theme_light()
#
#   # baseline
#   baseline <- sim_baseline_biv(object)
#   X1 <- baseline$X1[time_index,]
#   X2 <- baseline$X2[time_index,]
#   u2 <- copula::pobs(cbind(X1,X2))
#   kde.fit2 <- kdecop(u2)
#
#   contour_dat_P <- reshape2::melt(kde.fit2$estimate) %>%
#     mutate(u1 = kde.fit2$grid[Var1],
#            u2 = kde.fit2$grid[Var2]) %>%
#     select(u1,u2,value)
#
#   plot2 <- ggplot(contour_dat_P,aes(x=u1,y=u2,z=value)) +
#     geom_contour(aes(colour = after_stat(level)),
#                  breaks = c(0,0.5,1.5,1,2,3,4,5,10,15))+
#     scale_colour_gradient(low = "blue", high = "cyan")+ theme_light()
#
#   combined_plot <- plot1 +
#     scale_colour_gradient(low = "red", high = "mediumorchid1", name = "Q (level)") +
#     ggnewscale::new_scale_color() +
#     geom_contour(data=contour_dat_P,aes(colour = after_stat(level)),
#                  breaks = c(0,0.5,1,1.5,2,3,4,5,10,15),lty=2)+
#     scale_colour_gradient(low = "blue", high = "cyan", name = "P (level)")
#
#   return(list(plot1=plot1,plot2=plot2,combined_plot=combined_plot))
# }
#
# res1 <- plot_copula_dens_contour_v2(gamma_ind_ex,1)
# res1$combined_plot +
#   ggtitle("Copula density for terminal path values, jumps independent")
# ggpubr::ggarrange(res1$plot2+ggtitle("P-measure"),
#                   res1$plot1+ggtitle("Q-measure"))
#
# res2 <- plot_copula_dens_contour_v2(gamma_t_ex,1)
# res2$combined_plot +
#   ggtitle("Copula density for terminal path values, jumps independent")
# ggpubr::ggarrange(res2$plot2+ggtitle("P-measure"),
#                   res2$plot1+ggtitle("Q-measure"))

# res1 <- plot_copula_contour(gamma_ind_ex,1)
# res2 <- plot_copula_contour(gamma_t_ex,1)
#
#
# headings <- c('t' = "t copula", 'ind' = "Indepedent copula")
#
# rbind(dplyr::mutate(res1$data,copula="ind"),
#       dplyr::mutate(res2$data,copula="t")) %>%
#   ggplot2::ggplot(ggplot2::aes(x=x,y=y,z=z,colour=measure)) +
#   ggplot2::geom_contour(lwd=1.1) +
#   ggplot2::theme_bw(base_size = 14) +
#   ggplot2::facet_wrap(~copula,labeller = ggplot2::as_labeller(headings)) +
#   ggplot2::theme(legend.position = "bottom") +
#   ggplot2::xlab("X1") + ggplot2::ylab("X2")
#
#
#
# res1 <- plot_copula_contour(gamma_ind_ex,1)
# res2 <- plot_copula_contour(gamma_t_ex,1)
# res3 <- plot_copula_contour(gamma_counter_ex,1)
# res4 <- plot_copula_contour(gamma_comono_ex,1)
#
# headings <- c('t' = "t copula", 'ind' = "Indpedent copula",
#               "counter" = "Countermonotonic copula",
#               "comono" = "Comonotonic copula")
#
# rbind(dplyr::mutate(res1$data,copula="ind"),
#       dplyr::mutate(res2$data,copula="t"),
#       dplyr::mutate(res3$data,copula="counter"),
#       dplyr::mutate(res3$data,copula="comono")) %>%
#   ggplot2::ggplot(ggplot2::aes(x=x,y=y,z=z,colour=measure)) +
#   ggplot2::geom_contour(lwd=1) +
#   ggplot2::theme_minimal() +
#   ggplot2::facet_wrap(~copula,labeller = ggplot2::as_labeller(headings)) +
#   ggplot2::theme(legend.position = "bottom") +
#   ggplot2::xlab("X1") + ggplot2::ylab("X2")
#
#
# #
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
