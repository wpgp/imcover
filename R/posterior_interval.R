
#' Summarise posterior estimates of immunisation coverage
#'
#' @description Computes the Bayesian posterior summary intervals, sometimes
#'   called 'credible intervals' for the estimated immunisation coverage term
#'   ('mu').
#' @param X A fitted model object of type \code{icfit} or \code{iclist}.
#' @param object Character specifying the name of the object to summarise.
#'   Default is 'posterior'.
#' @param stat Character vector of summary statistics to apply.
#' @param probs Numeric vector of probabilities with values in [0,1] to be used
#'   for \code{stat} "quantile".
#' @return A \code{data.frame} with columns for 'country', 'time', and 'vaccine'
#'   labels and summary statistics, and as many rows as there are 'mu'
#'   observations of immunisation coverage.
#'
#' @name ic_coverage
#' @export
ic_coverage <- function(X,
                        object = 'posterior',
                        stat = c('mean', 'median', 'sd', 'quantile'),
                        probs = c(0.025, 0.25, 0.5, 0.75, 0.975)){
  UseMethod("ic_coverage")
}

#' @name ic_coverage
#' @export
ic_coverage.icfit <- function(X,
                              object = 'posterior',
                              stat = c('mean', 'median', 'sd', 'quantile'),
                              probs = c(0.025, 0.25, 0.5, 0.75, 0.975)){

  # match.arg(stat)
  stopifnot(all(probs >= 0 & probs <= 1))

  X <- X[[object]]

  psummary <- lapply(stat, FUN = function(st){
    if(st == 'quantile'){
      t(apply(X[,-c(1:3)], 1, st, probs))
    } else{
      apply(X[,-c(1:3)], 1, st)
    }
  })

  psummary <- do.call(cbind.data.frame, psummary)

  summarynm <- stat
  if('quantile' %in% stat){
    summarynm <- stat[-which(stat == 'quantile')]
    summarynm <- c(summarynm, paste0(probs * 100, '%'))
  }
  names(psummary) <- summarynm

  mu_hat <- cbind(X[, c(1:3)], psummary)
  return(mu_hat)
}


#' @name ic_coverage
#' @export
ic_coverage.iclist <- function(X,
                               object = 'posterior',
                               stat = c('mean', 'median', 'sd', 'quantile'),
                               probs = c(0.025, 0.25, 0.5, 0.75, 0.975)){
  out <- lapply(X, function(Xpost){
    ic_coverage(Xpost, object, stat, probs)
  })

  return(out)
}

