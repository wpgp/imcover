
#' Plot ic fit objects
#'
#' Provides methods to plot estimated time trends of vaccine coverage for
#' countries.
#' @param X Object of type \code{icfit} or \code{iclist}
#' @param observed Object of type \code{ic.df}
#' @param prob Numeric vector of length 2 for probabilities with values in [0,1]
#'   to be used to calculate the upper and lower credible intervals.
#' @param vaccine Optional. Character vector to subset which vaccines to plot
#' @param filter_yovi Logical. Should the estimates be filtered by year of
#'   vaccine introduction?
#' @param yovi Table containing years of vaccine introduction. WHO standard
#'   tables are used by default.
#' @return  A plot of the coverage estimates extracted from the fit, including
#'   credible intervals.
#'
#' @details If observed data are provided, they will be overlaid on the
#'   estimates.
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon geom_point
#' scale_x_continuous element_blank ggtitle theme_bw
#' @importFrom rlang .data
#'
#' @name ic_plot
#' @export
ic_plot <- function(X, observed,
                    prob = c(0.025, 0.975),
                    vaccine,
                    filter_yovi = TRUE, yovi){
  UseMethod("ic_plot")
}


#' @name ic_plot
#' @export
ic_plot.icfit <- function(X, observed,
                          prob = c(0.025, 0.975),
                          vaccine,
                          filter_yovi = TRUE, yovi){

  stopifnot(length(prob) == 2)

  # find years to plot
  tt <- ggplot2::cut_interval(list_times(X), length = 5)  # 5-year groups
  # clean labels
  lbls <- levels(tt)
  lbls <- paste(lbls, collapse = ',')
  lbls <- unique(as.numeric(strsplit(gsub("\\[|\\]|\\(", "", lbls), ',')[[1]]))

  # select observed data
  if(!missing(observed)){
    observed <- swap_names(observed)
    observed <- data.frame(observed)
  }

  # get the summary of the estimates ('mu')
  mu_hat <- posterior_interval(X, prob = prob)
  # find names to plot
  mu_hat[['lo']] <- mu_hat[[paste0(prob[1] * 100, '%')]]
  mu_hat[['hi']] <- mu_hat[[paste0(prob[2] * 100, '%')]]

  # filter by year of introduction
  if(filter_yovi){

  }

  # find vaccine subsets
  if(missing(vaccine)){
    vaccine <- list_vaccines(X)
  } else{
    vaccine <- intersect(vaccine, list_vaccines(X))
  }

  for(v in vaccine){
    plot <- ggplot2::ggplot() +
      ggplot2::geom_line(data = mu_hat[mu_hat$vaccine == v, ],
                         ggplot2::aes(x = .data$time, y = mean),
                         color = "black") +
      ggplot2::geom_ribbon(data = mu_hat[mu_hat$vaccine == v, ],
                           ggplot2::aes(x = .data$time, ymin = lo, ymax = hi),
                           alpha = 0.2, color = 'grey50')

    if(!missing(observed)){
      plot <- plot +
        ggplot2::geom_point(data = observed[observed$vaccine == v, ],
                            ggplot2::aes(x = .data$time, y = coverage,
                                         color = factor(source)))
    }

    plot <- plot +
      ggplot2::facet_wrap(. ~ country, scale = 'free', ncol = 4) +
      ggplot2::ggtitle(v) +
      ggplot2::scale_x_continuous(breaks = lbls,
                                  labels = lbls) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = 'bottom',
                     legend.title = ggplot2::element_blank())

    print(plot)
  }
}


#' @name ic_plot
#' @export
ic_plot.iclist <- function(X, observed,
                           prob = c(0.025, 0.975),
                           vaccine,
                           filter_yovi = TRUE, yovi){

  for(i in X){
    ic_plot(i, observed, prob, vaccine, filter_yovi, yovi)
  }
}
