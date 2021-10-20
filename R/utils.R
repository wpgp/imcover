#' List unique attributes
#'
#' Retrieve the unique list of attributes observed in an ic dataset.
#' @param X Object of class \code{ic.df}.
#' @name list_vaccines
#' @export
list_vaccines <- function(X){
  UseMethod("list_vaccines")
}

#' @name list_vaccines
#' @export
list_vaccines.ic.df <- function(X){
  return(sort(unique(X[[attr(X, "vaccine")]])))
}

#' @name list_vaccines
#' @export
list_vaccines.icfit <- function(X){
  return(sort(unique(X[['posterior']][['vaccine']])))
}


#' @rdname list_vaccines
#' @export
list_countries <- function(X){
  UseMethod("list_countries")
}

#' @rdname list_vaccines
#' @export
list_countries.ic.df <- function(X){
  return(sort(unique(X[[attr(X, "country")]])))
}

#' @rdname list_vaccines
#' @export
list_countries.icfit <- function(X){
  return(sort(unique(X[['posterior']][['country']])))
}


#' @rdname list_vaccines
#' @export
list_times <- function(X){
  UseMethod("list_times")
}

#' @rdname list_vaccines
#' @export
list_times.ic.df <- function(X){
  return(sort(unique(X[[attr(X, "time")]])))
}

#' @rdname list_vaccines
#' @export
list_times.icfit <- function(X){
  return(sort(unique(X[['posterior']][['time']])))
}

#' @rdname list_vaccines
#' @export
list_sources <- function(X){
  UseMethod("list_sources")
}

#' @rdname list_vaccines
#' @export
list_sources.ic.df <- function(X){
  return(sort(unique(X[[attr(X, "source")]])))
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


#' Convert attributes to column names
#' @keywords internal
swap_names <- function(X){
  stopifnot(is.ic_data(X))

  nn <- get_attr(X, ic_core())
  X <- data.frame(X)
  idx <- match(nn, names(X))
  names(X)[na.omit(idx)] <- names(nn)[which(!is.na(idx))]
  # convert back
  X <- ic_data(X,
               region = 'region', country = 'country',
               time = 'time', vaccine = 'vaccine',
               coverage = 'coverage', source = 'source')

  return(X)
}


#' Logit calculations
#' Helper function to calculate the logit and inverse-logit transforms.
#' @param x Numeric vector.
#' @return A numeric vector the same length as \code{x}.
#' @name logit
#' @export
logit <- function(x) stats::qlogis(x)


#' @rdname logit
#' @export
invlogit <- function(x) stats::plogis(x)

