
#' Summarise posterior estimates of immunisation coverage
#'
#' @description Computes the Bayesian posterior summary intervals, sometimes
#'   called 'credible intervals' for the estimated immunisation coverage term
#'   ('mu').
#' @param X A fitted model object of type \code{icfit} or \code{iclist}.
#' @param pars The parameters to calculate the intervals for.
#' @param stat Character vector of summary statistics to apply.
#' @param prob Numeric vector of probabilities with values in [0,1] to be used
#'   for \code{stat} "quantile".
#' @return A \code{data.frame} with columns for 'country', 'time', and 'vaccine'
#'   labels and summary statistics, and as many rows as there are 'mu'
#'   observations of immunisation coverage.
#'
#' @aliases posterior_interval
#' @importFrom rstantools posterior_interval
#' @export
#' @method posterior_interval icfit
#' @rdname posterior_interval.icfit
posterior_interval.icfit <- function(X,
                                     pars = 'mu',
                                     stat = c('mean', 'median', 'sd', 'quantile'),
                                     prob = c(0.025, 0.25, 0.5, 0.75, 0.975)){

  if(length(pars) > 1) pars <- pars[1L]

  if(pars != 'mu'){
    draws <- rstan::extract(X$fit, pars = pars)
    interval <- rstantools::posterior_interval(draws, prob)
    return(interval)
  }

  match.arg(stat)
  stopifnot(all(prob >= 0 & prob <= 1))

  X <- X[['posterior']]

  psummary <- lapply(stat, FUN = function(st){
    if(st == 'quantile'){
      t(apply(X[,-c(1:3)], 1, st, prob))
    } else{
      apply(X[,-c(1:3)], 1, st)
    }
  })

  psummary <- do.call(cbind.data.frame, psummary)

  summarynm <- stat
  if('quantile' %in% stat){
    summarynm <- stat[-which(stat == 'quantile')]
    summarynm <- c(summarynm, paste0(prob * 100, '%'))
  }
  names(psummary) <- summarynm

  mu_hat <- cbind(X[, c(1:3)], psummary)
  return(mu_hat)
}


#' @aliases posterior_interval
#' @importFrom rstantools posterior_interval
#' @export
#' @method posterior_interval iclist
#' @rdname posterior_interval.icfit
posterior_interval.iclist <- function(X,
                                      pars = 'mu',
                                      stat = c('mean', 'median', 'sd', 'quantile'),
                                      prob = c(0.025, 0.25, 0.5, 0.75, 0.975)){
  out <- lapply(X, function(Xpost){
    posterior_interval(Xpost, pars, stat, prob)
  })

  return(out)
}

