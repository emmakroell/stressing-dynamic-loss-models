#' plot a RPS object
#' choose to plot paths only, or compare to baseline model,
#' or compare to kappa_Q
#'
#' @param object RPS_model object
#' @param Npaths int, number of paths to plot
#' @param plot_type one of "compare_to_baseline" or "compare_to_kappa". If empty,
#' plots only stressed model
#' @param quantiles list containing lower and upper, specifies which quantiles
#' to highlight
#'
#' @return
#' @export
#'
path_plot <- function(object, Npaths, plot_type = "paths",
                      quantiles=list(lower=0.1,upper=0.9)) {

  X <- as.data.frame(object$paths)
  colnames(X) <- seq(1,dim(X)[2],1)
  time <- object$time_vec

  # select paths to plot randomly
  withr::with_seed(45,{
    indices <- sample(1:dim(X)[2],Npaths)
  })

  # reformat X for ggplot
  X_for_plot <- dplyr::as_tibble(cbind(X[,indices],time)) %>%
    tidyr::pivot_longer(cols = tidyselect::all_of(1:Npaths),
                        values_to = 'value',
                        names_to = 'number') %>%
    dplyr::mutate(number = as.factor(number), type = "stressed")

  # summarize quantiles of X for plot
  X_quantiles <- dplyr::as_tibble(cbind(X,time)) %>%
    tidyr::pivot_longer(cols = tidyselect::all_of(1:dim(X)[2]),
                        values_to = 'value',
                        names_to = 'number') %>%
    dplyr::group_by(time) %>%
    dplyr::summarise(lwr = stats::quantile(value,quantiles$lower,type=1,na.rm=T),
                     median = stats::median(value,na.rm=T),
                     upr = stats::quantile(value,quantiles$upper,type=1,na.rm=T)) %>%
    tidyr::pivot_longer(cols="lwr":"upr",values_to = "value",names_to ="quantile") %>%
    dplyr::mutate(quantile = as.factor(quantile),type = "stressed")

  if (plot_type == "paths"){
    plot <- ggplot2::ggplot() +
      ggplot2::geom_line(data=X_for_plot,
                         ggplot2::aes(time,value,colour=number),alpha=0.8) +
      ggplot2::geom_line(data=X_quantiles,
                         ggplot2::aes(x=time,y=value,group=quantile), size = 1) +
      ggplot2::theme_bw() + ggplot2::theme(legend.position = "none")

  } else if (plot_type == "compare_to_baseline") {
    X_baseline <- sim_baseline(object) %>% as.data.frame()
    colnames(X_baseline) <- seq(1,dim(X_baseline)[2])

    X_baseline_for_plot <- dplyr::as_tibble(cbind(X_baseline[,indices],time)) %>%
      tidyr::pivot_longer(cols = tidyselect::all_of(1:Npaths),
                          values_to = 'value',
                          names_to = 'number') %>%
      dplyr::mutate(number = as.factor(number), type = "original")

    # summarize quantiles of X for plot
    X_baseline_quantiles <- dplyr::as_tibble(cbind(X_baseline,time)) %>%
      tidyr::pivot_longer(cols = tidyselect::all_of(1:dim(X_baseline)[2]),
                   values_to = 'value', names_to = 'number') %>%
      dplyr::group_by(time) %>%
      dplyr::summarise(lwr = stats::quantile(value,quantiles$lower,type=1,na.rm=T),
                       median = stats::median(value,na.rm=T),
                       upr = stats::quantile(value,quantiles$upper,type=1,na.rm=T)) %>%
      tidyr::pivot_longer(cols="lwr":"upr",values_to = "value",names_to ="quantile") %>%
      dplyr::mutate(quantile = as.factor(quantile),type = "original")

    # headings for facet plot
    headings <- c('stressed'= 'stressed ~ X * -paths',
                  "original"= 'original ~ X * -paths')

    plot <- ggplot2::ggplot() +
      ggplot2::geom_line(data=rbind(X_for_plot,X_baseline_for_plot),
                         ggplot2::aes(time,value,colour=number),alpha=0.8) +
      ggplot2::geom_line(data=rbind(X_quantiles,X_baseline_quantiles),
                         ggplot2::aes(x=time,y=value,group=quantile), size = 1) +
      ggplot2::facet_wrap(~type,scales="free",
                          labeller = ggplot2::as_labeller(headings,ggplot2::label_parsed)) +
      ggplot2::theme_bw() + ggplot2::theme(legend.position = "none")


  } else if (plot_type == "compare_to_kappa") {
    kappa <- as.data.frame(object$kappa_Q)
    colnames(kappa) <- seq(1,dim(kappa)[2],1)

    kappa_for_plot <- dplyr::as_tibble(cbind(kappa[2:dim(kappa)[1],indices],
                                             time=time[2:length(time)])) %>%
      tidyr::pivot_longer(cols = tidyselect::all_of(1:Npaths),
                          values_to = 'value',
                          names_to = 'number') %>%
      dplyr::mutate(number = as.factor(number), type = "kappa")

    # headings for facet plot
    headings <- c('kappa'= 'stressed ~ kappa * -paths',
                  "stressed"= 'stressed ~ X * -paths')

    # horizontal lines for kappa^P and q
    lines <- dplyr::as_tibble(cbind(line=c(object$kappa,object$stress_parms$q),
                              type = c("kappa","stressed"))) %>%
             dplyr::mutate(line=as.numeric(line))

    plot <- rbind(X_for_plot,kappa_for_plot) %>%
      ggplot2::ggplot(ggplot2::aes(x=time,y=value,colour=number)) +
      ggplot2::geom_line(alpha=0.8) +
      ggplot2::facet_wrap(~type,scales="free",
                          labeller = ggplot2::as_labeller(headings,ggplot2::label_parsed)) +
      ggplot2::geom_hline(data = lines, ggplot2::aes(yintercept = line),lty=2) +
      ggplot2::theme_bw() + ggplot2::theme(legend.position = "none")

  } else{
    stop("Plot type not initialized")
  }
  return(plot)
}


