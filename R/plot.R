
#' Plot ic fit objects
#'
#' Provides methods to plot estimated time trends of vaccine coverage for
#' countries.
#' @param X Object of type \code{icfit} or \code{iclist}
#' @param ... Additional arguments passed to \code{ic_plot}.
#' @return Plot
#'
#' @seealso \code{\link[imcover]{ic_plot}}
#' @export
plot.icfit <- function(X, ...){
  ic_plot(X, ...)
}


#' Plot ic fit objects
#'
#' Provides methods to plot estimated time trends of vaccine coverage for
#' countries.
#' @param X Object of type \code{icfit} or \code{iclist}
#' @param ... Additional arguments passed to \code{ic_plot}.
#' @return Plot
#'
#' @seealso \code{\link[imcover]{ic_plot}}
#' @export
plot.iclist <- function(X, ...){
  ic_plot(X, ...)
}


#' Plot ic fit objects
#'
#' Provides methods to plot estimated time trends of vaccine coverage for
#' countries.
#' @param X Object of type \code{icfit} or \code{iclist}
#' @param probs Numeric vector of length 2 for probabilities with values in [0,1]
#'   to be used to calculate the upper and lower credible intervals.
#' @param vaccine Optional. Character vector to subset which vaccines to plot
#' @param observed Logical. Should the observed data used to fit \code{X} be
#'   overlaid on the plot? Default is \code{TRUE}.
#' @param prediction Logical. If predictions are found in the \code{icfit}
#'   object, should they be overlaid on the plot? Default is \code{TRUE}.
#' @param wuenic Logical. Should WUENIC vaccination coverage estimates be
#'   overlaid on the plot? Default is \code{FALSE}.
#' @param filter_yovi Logical. Should the estimates be filtered by year of
#'   vaccine introduction? Default is \code{FALSE}.
#' @param ncol When plotting multiple countries, the number of plots per page.
#' @param interactive Logical. Should the output be interactive? Default is
#'   \code{FALSE}.
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
ic_plot <- function(X,
                    probs = c(0.025, 0.975),
                    vaccine,
                    observed = TRUE,
                    prediction = TRUE,
                    wuenic = FALSE,
                    filter_yovi = FALSE,
                    ncol = 4,
                    interactive = FALSE){
  UseMethod("ic_plot")
}


#' @name ic_plot
#' @export
ic_plot.icfit <- function(X,
                          probs = c(0.025, 0.975),
                          vaccine,
                          observed = TRUE,
                          prediction = TRUE,
                          wuenic = FALSE,
                          filter_yovi = FALSE,
                          ncol = 4,
                          interactive = FALSE){

  stopifnot(length(probs) == 2)

  # select observed data
  if(observed){
    obsdat <- swap_names(X$data)
    obsdat <- data.frame(obsdat)
    obsdat$source <- factor(obsdat$source)
    if(is.null(obsdat)) observed <- FALSE
  }

  # get the summary of the estimates ('mu')
  X[['posterior']] <- ic_coverage(X, "posterior", probs = probs)

  # filter by year of introduction
  if(filter_yovi){
    X <- filter_yovi(X, na.rm = TRUE)
  }

  # extract posterior values
  mu_hat <- X$posterior

  # set-up for predictions
  if(prediction){
    if(!is.null(X$prediction)){
      pred <- ic_coverage(X, "prediction", probs = probs)
      t0 <- min(pred$time)

      mu_hat <- rbind(mu_hat, pred)
    } else{
      prediction <- FALSE
    }
  }

  # set-up for WUENIC data
  if(wuenic){
    wd <- download_wuenic(return_ic = FALSE, quiet = TRUE)
    # subset cols
    wd <- wd[, c("ISOCountryCode", "Year", "Vaccine", "WUENIC")]
    wd$ISOCountryCode <- toupper(wd$ISOCountryCode)
    wd$Vaccine <- toupper(wd$Vaccine)
    names(wd) <- c('country', 'time', 'vaccine','coverage')
    wd$source <- "WUENIC"

    # filter
    wd <- wd[wd$country %in% list_countries(X), ]
    wd <- wd[wd$time %in% list_times(X), ]
    wd <- wd[wd$vaccine %in% list_vaccines(X), ]

    # check
    if(nrow(wd) == 0L){
      wuenic <- FALSE
    }
  }

  # find years to plot
  tt <- ggplot2::cut_interval(mu_hat$time, length = 5)  # 5-year groups
  # clean labels
  lbls <- levels(tt)
  lbls <- paste(lbls, collapse = ',')
  lbls <- unique(as.numeric(strsplit(gsub("\\[|\\]|\\(", "", lbls), ',')[[1]]))

  # find names to plot
  mu_hat[['lo']] <- mu_hat[[paste0(probs[1] * 100, '%')]]
  mu_hat[['hi']] <- mu_hat[[paste0(probs[2] * 100, '%')]]

  # find vaccine subsets
  if(missing(vaccine)){
    vaccine <- list_vaccines(X)
  } else{
    vaccine <- intersect(vaccine, list_vaccines(X))
  }

  for(v in vaccine){
    plotobj <- ggplot2::ggplot() +
      ggplot2::geom_line(data = mu_hat[mu_hat$vaccine == v, ],
                         ggplot2::aes(x = .data$time,
                                      y = .data$mean),
                         color = "black") +
      ggplot2::geom_ribbon(data = mu_hat[mu_hat$vaccine == v, ],
                           ggplot2::aes(x = .data$time,
                                        ymin = .data$lo, ymax = .data$hi),
                           alpha = 0.2, color = 'grey50')

    if(prediction){
      plotobj <- plotobj +
        ggplot2::geom_vline(xintercept = t0, lty = 'dashed', color = 'blue')
    }

    if(observed){
      plotobj <- plotobj +
        ggplot2::geom_point(data = obsdat[obsdat$vaccine == v & obsdat$country %in% mu_hat$country, ],
                            ggplot2::aes(x = .data$time,
                                         y = .data$coverage,
                                         color = source))
    }

    if(wuenic){
      plotobj <- plotobj +
        ggplot2::geom_point(data = wd[wd$vaccine == v & wd$country %in% mu_hat$country, ],
                            ggplot2::aes(x = .data$time,
                                         y = .data$coverage,
                                         color = source),
                            shape = 17)
    }

    plotobj <- plotobj +
      ggplot2::facet_wrap(. ~ country, scale = 'free', ncol = ncol) +
      ggplot2::ggtitle(v) +
      ggplot2::scale_x_continuous(breaks = lbls,
                                  labels = lbls) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = 'bottom',
                     legend.title = ggplot2::element_blank())

    # output
    if(interactive){
      print(plotly::ggplotly(plotobj))
    } else{
      print(plotobj)
    }
  }
}


#' @name ic_plot
#' @export
ic_plot.iclist <- function(X,
                           probs = c(0.025, 0.975),
                           vaccine,
                           observed = TRUE,
                           prediction = TRUE,
                           filter_yovi = FALSE,
                           ncol = 4,
                           interactive = FALSE){

  for(i in X){
    ic_plot(i, probs, vaccine, observed, prediction,
            filter_yovi, ncol, interactive)
  }
}
