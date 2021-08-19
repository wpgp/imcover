#' List vaccines
#'
#' Retrieve the unique list of vaccines observed in an ic dataset.
#' @param X Object of class \code{ic.df}.
#' @name list_vaccines
#' @export
list_vaccines <- function(X){
  if(!is.ic_data(X)){
    stop("Please supply valid 'ic' data")
  }
  return(unique(X[[attr(X, "vaccine")]]))
}


#' List core attributes for ic data
#' Retrieve the names of the core attributes that can be expected in
#' \code{ic.df} data.
#' @param survey Include survey-specific attributes or only \code{ic.df}
#'   attributes? Default is \code{TRUE}.
#' @name ic_core
#' @export
ic_core <- function(survey = TRUE){
  stopifnot(is.logical(survey))
  if(survey){
    return(c("region", "country", "time", "vaccine", "coverage", "source", "dose",
             "population", "source", "survey", "evidence", "validity", "sample"))
  } else{
    return(c("region", "country", "time", "vaccine", "coverage",
             "source", "dose", "population"))
  }
}


#' Get the attributes of ic data
#' Retrieve the value(s) of attributes from an \code{ic.df} dataset. While
#' intended to be used with the core ic attributes, it will work with any named
#' attribute.
#' @param X An object of class \code{ic.df}.
#' @param attrs A vector of character names of attributes to be retrieved.
#'   Default are all core attributes.
#' @param unlist Should the function return a named list of attributes and
#'   values or a named vector? Default is \code{TRUE} to return a vector.
#' @name get_attr
#' @export
get_attr <- function(X, attrs = ic_core(), unlist = TRUE){
  if(unlist){
    return(unlist(attributes(X)[which(names(attributes(X)) %in% attrs)]))
  } else{
    return(attributes(X)[which(names(attributes(X)) %in% attrs)])
  }
}


#' Get year of introduction for vaccines
#' Retrieve a list of the year when a vaccine was introduced to the entire
#' country.
#' @param vaccine Vector of vaccine codes as characters. Default is \code{NULL}
#'   which returns all vaccines in the dataset.
#' @details Data source of vaccine schedules:
#'   \link{https://www.who.int/immunization/monitoring_surveillance/data/year_vaccine_introduction.xlsx}
#'
#' @name get_yovi
#' @export
get_yovi <- function(vaccine = NULL){
  if(!is.null(vaccine)){
    return(imcover:::yovi[vaccine %in% yovi$vaccine, ])
  } else{
    return(imcover:::yovi)
  }
}


#' Filter ic data records based on year of vaccine introduction
#' Exclude ic data records for years when a vaccine had not been fully
#' introduced in a country.
#' @param X An object of class \code{ic.df}.
#' @param vaccine Vector of vaccine codes as characters. Default is \code{NULL}
#'   which includes all vaccines in object \code{X}.
#' @param yovi Table of the year of vaccine introduction. See details. Default
#'   value comes from \code{imcover::get_yovi()}.
#' @param na.rm Should records which do not have any data for year of
#'   introduction be removed? Default is \code{FALSE} to keep \code{NA} records
#'   in \code{X}.
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
  if(!is.ic_data(X)){
    stop("Please provide valid 'ic' data.")
  }

  if(!inherits(yovi, "data.frame")){
    stop("Please provide valid 'yovi' data.")
  }

  if(any(!c('ISO3code', 'vaccine', 'year_introduced') %in% names(yovi))){
    stop("Non-matching column names in YOVI table.")
  }

  if(!is.null(vaccine)){
    yovi <- yovi[yovi$vaccine %in% vaccine, ]
  } else{
    vaccine <- unique(yovi$vaccine)
  }
  if(nrow(yovi) < 1L) stop("No vaccine records found in YOVI table.")

  X <- merge(X, yovi[, c("ISO3code", "vaccine", "year_introduced")],
             by.x = get_attr(X, attrs = c("country", "vaccine")),
             by.y = c("ISO3code", "vaccine"),
             suffixes = c("", ".yovi"),
             all.x = TRUE, sort = FALSE)
  nm <- c("year_introduced.yovi", "year_introduced")
  nm <- nm[which(nm %in% names(X))]

  # filter
  if(na.rm){
    del <- which((is.na(X[[nm]]) & X[[get_attr(X, "vaccine")]] %in% vaccine) | X[[get_attr(X, "time")]] <= X[[nm]])
  } else{
    del <- which(!is.na(X[[nm]]) & X[[get_attr(X, "time")]] <= X[[nm]])
  }
  X[[nm]] <- NULL

  return(X[-del, ])
}


#' Get regional groupings for country codes
#' Retrieve a list of the region based on the ISO3c alpha codes.
#' @param country Vector of country ISO3 codes as characters. Default is \code{NULL}
#'   which returns all regions in the dataset.
#' @param type The type of region grouping to return. The default \code{'who'}
#'   returns World Health Organisation codes. The other option is 'm49' for UN
#'   Stats Division region groupings.
#' @details Data source of M49 regions:
#'   \link{https://unstats.un.org/unsd/methodology/m49/overview/}
#' @return A character vector of region codes of length equal
#'   \code{length(country)} if \code{country} is specified or a \code{data.frame} of
#'   the region table if \code{country} is \code{NULL}.
#'
#' @name get_region
#' @export
get_region <- function(country = NULL, type = 'who'){

  if(type == 'm49'){
    region <- imcover:::m49region
  } else if(type == 'who'){
    region <- imcover:::whoregion
  } else{
    stop('Cannot find specified region table.')
  }

  if(!is.null(country)){
    return(merge(data.frame('ISO3code' = country),
                 region,
                 by = 'ISO3code', all.x = T)[,'region'])
  } else{
    return(region)
  }
}


#' @export
#' @param x Numeric vector.
#' @return A numeric vector the same length as \code{x}.
logit <- function(x) stats::qlogis(x)


#' @rdname logit
#' @export
invlogit <- function(x) stats::plogis(x)

