
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
#'   affect the posterior estimates and prediction (if present). When a
#'   selection parameter is omitted, then all possible values will be selected.
#'   Selection parameters which are not present in \code{X} will be ignored.
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
  pred <- X[['prediction']]

  if(!missing(country)){
    post <- post[post$country %in% country, ]

    if(!is.null(pred)){
      pred <- pred[pred$country %in% country, ]
    }

    if(nrow(post) == 0) return(NULL)
  }

  if(!missing(time)){
    post <- post[post$time %in% time, ]

    if(!is.null(pred)){
      pred <- pred[pred$time %in% time, ]
    }

    if(nrow(post) == 0) return(NULL)
  }

  if(!missing(vaccine)){
    post <- post[post$vaccine %in% vaccine, ]

    if(!is.null(pred)){
      pred <- pred[pred$vaccine %in% vaccine, ]
    }

    if(nrow(post) == 0) return(NULL)
  }

  out <- list('fit' = X$fit,
              'posterior' = post,
              'data' = X$data,
              'labels' = X$labels,
              'numerator' = X$numerator,
              'denominator' = X$denominator,
              'prediction' = pred)

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


#' Filter ic data records based on year of vaccine introduction
#' Exclude ic data records for years when a vaccine had not been fully
#' introduced in a country.
#' @param X An object of class \code{ic.df} or \code{icfit}.
#' @param vaccine Vector of vaccine codes as characters. Default is \code{NULL}
#'   which includes all vaccines in object \code{X}.
#' @param yovi Table of the year of vaccine introduction. See details. Default
#'   value comes from \code{imcover::get_yovi()}.
#' @param na.rm Should records which do not have any data for year of
#'   introduction be removed? Default is \code{FALSE} to keep \code{NA} records
#'   in \code{X}.
#' @return An object with type of \code{X} with records filtered by the year of
#'   vaccine introductions.
#' @details The \code{yovi} argument can take a user-defined table of values.
#'   This must be a \code{data.frame} object with column names 'ISO3code',
#'   'vaccine' and 'year_introduced'.
#' @seealso \code{\link[imcover]{get_yovi}}
#' @name filter_yovi
#' @export
filter_yovi <- function(X,
                        vaccine = NULL,
                        yovi = imcover::get_yovi(),
                        na.rm = FALSE){
  UseMethod("filter_yovi")
}


#' @name filter_yovi
#' @export
filter_yovi.ic.df <- function(X,
                              vaccine = NULL,
                              yovi = imcover::get_yovi(),
                              na.rm = FALSE){

  if(!inherits(yovi, "data.frame")){
    stop("Please provide valid 'yovi' data.")
  }

  if(any(!c('ISO3Code', 'vaccine', 'year_introduced') %in% names(yovi))){
    stop("Non-matching column names in YOVI table.")
  }

  if(!is.null(vaccine)){
    yovi <- yovi[yovi$vaccine %in% vaccine, ]
  } else{
    vaccine <- unique(yovi$vaccine)
  }
  if(nrow(yovi) < 1L) stop("No vaccine records found in YOVI table.")

  X <- merge(X, yovi[, c("ISO3Code", "vaccine", "year_introduced")],
             by.x = get_attr(X, attrs = c("country", "vaccine")),
             by.y = c("ISO3Code", "vaccine"),
             suffixes = c("", ".yovi"),
             all.x = TRUE, sort = FALSE)
  nm <- c("year_introduced.yovi", "year_introduced")
  nm <- nm[which(nm %in% names(X))]

  # filter
  if(na.rm){
    # del <- which((is.na(X[[nm]]) & X[[get_attr(X, "vaccine")]] %in% vaccine) | X[[get_attr(X, "time")]] <= X[[nm]])
    del <- ifelse((is.na(X[[nm]]) & X[[get_attr(X, "vaccine")]] %in% vaccine) | X[[get_attr(X, "time")]] <= X[[nm]], TRUE, FALSE)
  } else{
    # del <- which(!is.na(X[[nm]]) & X[[get_attr(X, "time")]] <= X[[nm]])
    del <- ifelse(!is.na(X[[nm]]) & X[[get_attr(X, "time")]] <= X[[nm]], TRUE, FALSE)
  }
  X[[nm]] <- NULL

  return(X[which(!del), ])
}


#' @name filter_yovi
#' @export
filter_yovi.icfit <- function(X,
                              vaccine = NULL,
                              yovi = imcover::get_yovi(),
                              na.rm = FALSE){

  if(!inherits(yovi, "data.frame")){
    stop("Please provide valid 'yovi' data.")
  }

  if(any(!c('ISO3Code', 'vaccine', 'year_introduced') %in% names(yovi))){
    stop("Non-matching column names in YOVI table.")
  }

  if(!is.null(vaccine)){
    yovi <- yovi[yovi$vaccine %in% vaccine, ]
  } else{
    vaccine <- unique(yovi$vaccine)
  }
  if(nrow(yovi) < 1L) stop("No vaccine records found in YOVI table.")

  # extract posterior data
  dat <- X$posterior

  dat <- merge(dat, yovi[, c("ISO3Code", "vaccine", "year_introduced")],
               by.x = c("country", "vaccine"),
               by.y = c("ISO3Code", "vaccine"),
               suffixes = c("", ".yovi"),
               all.x = TRUE,
               sort = FALSE)
  nm <- c("year_introduced.yovi", "year_introduced")
  nm <- nm[which(nm %in% names(dat))]

  # filter
  if(na.rm){
    # del <- which((is.na(dat[[nm]]) & dat[["vaccine"]] %in% vaccine) | dat[["time"]] <= dat[[nm]])
    del <- ifelse((is.na(dat[[nm]]) & dat[['vaccine']] %in% vaccine) | dat[['time']] <= dat[[nm]], TRUE, FALSE)
  } else{
    # del <- which(!is.na(dat[[nm]]) & dat[["time"]] <= dat[[nm]])
    del <- ifelse(!is.na(dat[[nm]]) & dat[['time']] <= dat[[nm]], TRUE, FALSE)
  }
  dat[[nm]] <- NULL

  X[['posterior']] <- dat[which(!del), ]

  return(X)
}


