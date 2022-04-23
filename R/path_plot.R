#' plot a RPS object
#' creates three-panel plot of paths under stressed model, paths under
#'  baseline model, and kappa_Q
#'
#' @param object RPS_model object
#' @param Npaths int, number of paths to plot
#' @param quantiles list containing lower and upper, specifies which quantiles
#' to highlight
#' @param indices vector of indices to plot of length Npaths,
#' or ="random" to choose randomly
#'
#' @return
#' @export
#'
path_plot_biv <- function(object, Npaths, quantiles=list(lower=0.1,upper=0.9),
                      indices = "random") {
  # recover time
  time <- object$time_vec

  # reshape X1
  X <- as.data.frame(object$paths$X1)
  colnames(X) <- seq(1,dim(X)[2],1)

  # select paths to plot randomly
  if (indices[1] == "random"){
    withr::with_seed(450,{
      indices <- sample(1:dim(X)[2],Npaths)
    })}

  # reformat X1 for ggplot
  X_for_plot <- dplyr::as_tibble(cbind(X[,indices],time)) %>%
    tidyr::pivot_longer(cols = tidyselect::all_of(1:Npaths),
                        values_to = 'value',
                        names_to = 'number') %>%
    dplyr::mutate(number = as.factor(number), type = "X1")

  # summarize quantiles of X1 for plot
  X_quantiles <- dplyr::as_tibble(cbind(X,time)) %>%
    tidyr::pivot_longer(cols = tidyselect::all_of(1:dim(X)[2]),
                        values_to = 'value',
                        names_to = 'number') %>%
    dplyr::group_by(time) %>%
    dplyr::summarise(lwr = stats::quantile(value,quantiles$lower,type=1,na.rm=T),
                     median = stats::median(value,na.rm=T),
                     upr = stats::quantile(value,quantiles$upper,type=1,na.rm=T)) %>%
    tidyr::pivot_longer(cols="lwr":"upr",values_to = "value",names_to ="quantile") %>%
    dplyr::mutate(quantile = as.factor(quantile),type = "X1")


  # X2
  X2 <- as.data.frame(object$paths$X2)
  colnames(X2) <- seq(1,dim(X2)[2],1)

  X2_for_plot <- dplyr::as_tibble(cbind(X2[,indices],time)) %>%
    tidyr::pivot_longer(cols = tidyselect::all_of(1:Npaths),
                        values_to = 'value',
                        names_to = 'number') %>%
    dplyr::mutate(number = as.factor(number), type = "X2")

  # summarize quantiles of X2 for plot
  X2_quantiles <- dplyr::as_tibble(cbind(X2,time)) %>%
    tidyr::pivot_longer(cols = tidyselect::all_of(1:dim(X2)[2]),
                        values_to = 'value', names_to = 'number') %>%
    dplyr::group_by(time) %>%
    dplyr::summarise(lwr = stats::quantile(value,quantiles$lower,type=1,na.rm=T),
                     median = stats::median(value,na.rm=T),
                     upr = stats::quantile(value,quantiles$upper,type=1,na.rm=T)) %>%
    tidyr::pivot_longer(cols="lwr":"upr",values_to = "value",names_to ="quantile") %>%
    dplyr::mutate(quantile = as.factor(quantile),type = "X2")

  # kappa
  kappa <- as.data.frame(object$kappa_Q)
  colnames(kappa) <- seq(1,dim(kappa)[2],1)

  kappa_for_plot <- dplyr::as_tibble(cbind(kappa[2:dim(kappa)[1],indices],
                                           time=time[2:length(time)])) %>%
  tidyr::pivot_longer(cols = tidyselect::all_of(1:Npaths),
                      values_to = 'value',
                      names_to = 'number') %>%
    dplyr::mutate(number = as.factor(number), type = "kappa")

  # order of panels:
  panel_order <- c("X1","kappa","X2")


  # horizontal lines for kappa^P and q
  lines <- dplyr::as_tibble(cbind(line=c(object$stress_parms$q,object$kappa,
                                         object$stress_parms$q),
                                  type = panel_order)) %>%
    dplyr::mutate(line=as.numeric(line),
                  type = factor(type,levels = panel_order))

  # headings for facet plot
  headings <- c('X1'= 'X[1] ~ paths',
                'X2'= 'X[2] ~ paths',
                'kappa'= 'kappa ~ paths')

  # enforce once scale
  max_X <- max(max(X_for_plot$value),max(X2_for_plot$value))
  blank_data <- rbind(dplyr::as_tibble(cbind(time=time,value=max_X,type="X1")),
                      dplyr::as_tibble(cbind(time=time,value=max_X,type="X2"))) %>%
    dplyr::mutate(time=as.numeric(time),
                  value=as.numeric(value),
                  type = factor(type,levels = panel_order))

  # make plot
  rbind(X_for_plot,X2_for_plot,kappa_for_plot) %>%
    dplyr::mutate(type = factor(type,levels = panel_order)) %>%
    ggplot2::ggplot() +
    ggplot2::geom_line(alpha=0.8,ggplot2::aes(x=time,y=value,colour=number),size=1) +
    ggplot2::geom_line(data=dplyr::mutate(rbind(X_quantiles,X2_quantiles),
                                   type = factor(type,levels = panel_order)),
                       ggplot2::aes(x=time,y=value,group=quantile), size = 1.1) +
    ggplot2::geom_hline(data = lines, ggplot2::aes(yintercept = line),lty=2) +
    ggplot2::facet_wrap(~type, scales="free",
                        labeller = ggplot2::as_labeller(headings,ggplot2::label_parsed)) +
    ggplot2::geom_blank(data=blank_data,ggplot2::aes(x=time,y=value)) +
    ggplot2::theme_bw(base_size=14) + ggplot2::theme(legend.position = "none") +
    ggplot2::ylab('')

}
