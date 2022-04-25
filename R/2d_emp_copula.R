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
#' @param type string, one of "copula", "mixture"
#'
#' @return
#' @export
#'
#' @examples
plot_copula <- function(object,time,type="copula"){
  time_index <- match(time,object$time_vec)

  # stressed
  X1 <- object$paths$X1[time_index,]
  X2 <- object$paths$X2[time_index,]
  stressed <- emp_copula(X1,X2) %>%
    dplyr::mutate(type = "stressed")

  # baseline
  if (type == "copula") baseline <- sim_baseline_biv(object)
  else if (type == "mixture") baseline <- sim_baseline_mixture(object)

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
}


#' Plot empirical copula contour plot under P and Q
#'
#' @param object RPS_model object
#' @param time time at which to plot
#' @param type string, one of "copula", "mixture"
#'
#' @return
#' @export
#'
#' @examples
plot_copula_contour <- function(object,time,type="copula",show_true_ind=FALSE){

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
  if (type == "copula") baseline_sim <- sim_baseline_biv(object)
  else if (type == "mixture") baseline_sim <- sim_baseline_mixture(object)

  X1 <- baseline_sim$X1[time_index,]
  X2 <- baseline_sim$X2[time_index,]
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


  if (show_true_ind == TRUE){
    xy <- NULL
    for (constant in seq(0.1,0.9,by=0.1)){
      xy <- rbind(xy, dplyr::as_tibble(list(x = seq(0,1,0.01),
                                            y = sapply(seq(0,1,0.01),
                                                       function(x) min(constant/x,1)),
                                            constant = constant)))
    }
    xy <- dplyr::mutate(xy,dplyr::across(3,as.numeric),
                        y = ifelse(y==1,NA,y))
  }


  # combined plot
  data <- rbind(res_P,res_Q)
  plot <- data %>%
    ggplot2::ggplot() +
    ggplot2::geom_contour(ggplot2::aes(x,y,z=z,colour=measure), lwd=1) +
    ggplot2::theme_minimal() +
    # ggplot2::geom_line(data = xy, ggplot2::aes(x,y,group=constant), lty = 2)  +
    ggplot2::scale_colour_manual(values=c("#00BFC4","#F8766D"))

  if (show_true_ind == TRUE) {
    plot <- plot +
      ggplot2::geom_line(data = xy, ggplot2::aes(x,y,group=constant), lty = 2)
  }

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

  if (object$model_type != "bivariate_copula"){
    stop(paste("Density contour plot method not defined for model type",
               object$model_type))
  }

  time_index <- match(time,object$time_vec)

  # stressed
  X1 <- object$paths$X1[time_index,]
  X2 <- object$paths$X2[time_index,]
  u1 <- copula::pobs(cbind(X1,X2))
  kde.fit1 <- kdecopula::kdecop(u1)
  plot(kde.fit1)
  graphics::par(mar=c(5,5,4,1)+.1)
  graphics::contour(kde.fit1,margins="unif",levels = c(0.1,0.25,0.5,1,1.5,2,3,4),
                    col="#F8766D",lwd=3,cex.lab=1.5,cex.axis=1.5,
                    labcex=1.1,vfont=c("sans serif", "bold"),
                    xlab=expression(X[1]),ylab=expression(X[2]))

  # baseline
  baseline <- sim_baseline_biv(object)
  X1 <- baseline$X1[time_index,]
  X2 <- baseline$X2[time_index,]
  u2 <- copula::pobs(cbind(X1,X2))
  kde.fit2 <- kdecopula::kdecop(u2)
  plot(kde.fit2)
  graphics::contour(kde.fit2,margins="unif",levels = c(0.1,0.25,0.5,1,1.5,2,3,4),
                    col="#00BFC4",add=TRUE,lwd=3,
                    labcex=1.1,vfont=c("sans serif", "bold"))
  # add legend
  graphics::legend("bottomright", legend = c("Q", "P"),
         col = c("#F8766D","#00BFC4"),lty=1:1,lwd=c(2,2),cex=1.1)
  graphics::par(mar=c(5,4,4,2)+.1)

}


