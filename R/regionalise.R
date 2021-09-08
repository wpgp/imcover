
#' Regional immunisation coverage
#'
#' Calculates population-weighted regional estimates from \code{icfit} objects.
#' @param X Object of \code{icfit} or \code{iclist} from model fitting
#' @param denom Population denominator data. Default is downloaded from WHO.
#' @param stat Character vector of summary statistics to apply.
#' @param probs Numeric vector of probabilities with values in [0,1] to be used
#'   for \code{stat} "quantile".
#' @param filter_yovi Logical. Should the estimates be filtered by year of
#'   vaccine introduction?
#' @return \code{data.frame} object with regional summaries
#' @name ic_regional
#' @export
ic_regional <- function(X,
                        denom,
                        stat = c('mean', 'median', 'sd', 'quantile'),
                        probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
                        filter_yovi = TRUE){
  UseMethod("ic_regional")
}


#' @name ic_regional
#' @export
ic_regional.icfit <- function(X,
                              denom,
                              stat = c('mean', 'median', 'sd', 'quantile'),
                              probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
                              filter_yovi = TRUE){

  if(missing(denom)){
    denom <- download_denom(use_cache = TRUE)
  } else{
    #check denom
  }

  if(filter_yovi){
    X <- filter_yovi(X, na.rm = TRUE)
  }

  # extract samples
  mu <- X[['posterior']]
  # inner join to select
  mu <- merge(mu,
              denom[, c('country', 'time', 'vaccine', 'population')],
              by = c('country', 'time', 'vaccine'))

  # split - apply - stack
  splits <- split(mu, mu[, c("time", "vaccine")])

  dat <- lapply(splits, function(sdat){
    totpop <- sum(sdat$population)
    sdat <- sdat$population * sdat[, 4:(ncol(sdat)-1)]
    return(colSums(sdat) / totpop)
  })

  dat <- do.call(rbind, dat)
  nmdat <- do.call(rbind.data.frame, strsplit(rownames(dat), ".", fixed = T))
  names(nmdat) <- c('year', 'vaccine')
  # clean up missings
  dat[is.nan(dat)] <- NA

  # recombine
  dat <- cbind(nmdat, dat)
  rownames(dat) <- NULL

  # summarise
  psummary <- lapply(stat, FUN = function(st){
    if(st == 'quantile'){
      t(apply(dat[,-c(1:2)], 1, st, probs, na.rm = TRUE))
    } else{
      apply(dat[,-c(1:2)], 1, st, na.rm = TRUE)
    }
  })

  psummary <- do.call(cbind.data.frame, psummary)

  summarynm <- stat
  if('quantile' %in% stat){
    summarynm <- stat[-which(stat == 'quantile')]
    summarynm <- c(summarynm, paste0(probs * 100, '%'))
  }
  names(psummary) <- summarynm

  dat <- cbind(dat[, c(1:2)], psummary)
  return(dat)
}


#' @name ic_regional
#' @export
ic_regional.iclist <- function(X,
                               denom,
                               stat = c('mean', 'median', 'sd', 'quantile'),
                               probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
                               filter_yovi = TRUE){

  out <- lapply(X, function(i){ ic_regional(i, stat = stat, probs = probs, filter_yovi = filter_yovi) })

  return(out)
}

