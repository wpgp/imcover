
#' Filter immunisation coverage records
#'
#' Select records from 'ic' objects based on region group, country, time, or
#' vaccine.
#' @param X object of type \code{ic.df}, \code{icfit}, or \code{iclist}
#' @param region Character vector of regional groups to select. Not valid for
#'   \code{icfit}.
#' @param country Character vector of country identifiers to select
#' @param time Vector of time points to select
#' @param vaccine Character vector of vaccine records to select, given by
#'   abbreviations (e.g. 'DTP1')
#' @return A new object with type matching \code{X} with the selected records.
#' @details For \code{icfit} and \code{iclist} objects, filtering will only
#'   affect the posterior estimates. When a selection parameter is omitted, then
#'   all possible values will be selected. Selection parameters which are not
#'   present in \code{X} will be ignored.
#'
#' @name ic_filter
#' @export
ic_filter <- function(X, region, country, time, vaccine) {
  UseMethod("ic_filter")
}


#' @name ic_filter
#' @export
ic_filter.ic.df <- function(X, region, country, time, vaccine){
  if(!missing(region)){
    X <- X[X[[get_attr(X, 'region')]] %in% region, ]
    if(nrow(X) == 0){ stop("No records found.") }
  }

  if(!missing(country)){
    X <- X[X[[get_attr(X, 'country')]] %in% country, ]
    if(nrow(X) == 0){ stop("No records found.") }
  }

  if(!missing(time)){
    X <- X[X[[get_attr(X, 'time')]] %in% time, ]
    if(nrow(X) == 0){ stop("No records found.") }
  }

  if(!missing(vaccine)){
    X <- X[X[[get_attr(X, 'vaccine')]] %in% vaccine, ]
    if(nrow(X) == 0){ stop("No records found.") }
  }

  return(X)
}


#' @name ic_filter
#' @export
ic_filter.icfit <- function(X, country, time, vaccine){
  post <- X[['posterior']]

  if(!missing(country)){
    post <- post[post$country %in% country, ]
    if(nrow(post) == 0) return(NULL)
  }

  if(!missing(time)){
    post <- post[post$time %in% time, ]
    if(nrow(post) == 0) return(NULL)
  }

  if(!missing(vaccine)){
    post <- post[post$vaccine %in% vaccine, ]
    if(nrow(post) == 0) return(NULL)
  }

  out <- list('fit' = X$fit, 'posterior' = post, 'data' = X$vax_data)
  class(out) <- list('icfit', class(out))

  return(out)
}


#' @name ic_filter
#' @export
ic_filter.iclist <- function(X, region, country, time, vaccine){
  if(!missing(region)){
    X <- X[names(X) %in% region]
    if(length(X) < 1L) return(NULL)
    if(is.list(X)) class(X) <- list('iclist', class(X))
    if(length(X) == 1L) X <- X[[1L]]
  }

  if(any(!missing(country), !missing(time), !missing(vaccine))){
    if(inherits(X, 'iclist')){
      out <- lapply(X, ic_filter, country = country, time = time, vaccine = vaccine)
      names(out) <- names(X)
      out <- out[lengths(out) != 0]  # drop NULLs
      class(out) <- list("iclist", class(out))
    } else{
      out <- ic_filter(X, country = country, time = time, vaccine = vaccine)
    }

    return(out)
  }

  return(X)
}

