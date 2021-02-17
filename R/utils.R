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
    return(c("group", "time", "vaccine", "coverage", "dose", "population",
             "survey", "evidence", "validity", "sample"))
  } else{
    return(c("group", "time", "vaccine", "coverage", "dose", "population"))
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


# ic_filter <- function(X, yovi = 'default'){
#   if(!is.ic_data(X)){
#     stop("Please provide valid 'ic' data.")
#   }
#
#   if(missing(yovi)){
#     stop("Please provide a data.frame of vaccines, or select 'default'.")
#   }
#
#   if(is.character(yovi)){
#     if(yovi == "default"){
#       yovi <- get_yovi()
#     }
#   } else{
#     stopifnot(inherits(yovi, "data.frame"))
#   }
#
#   # filter
#   return(X)
# }

#' Get year of introduction for vaccines
#' Retrieve a list of the year when a vaccine was introduced to the entire
#' country.
#' @param vaccine Vector of vaccine codes as characters. Default is \code{NULL}
#'   which returns all vaccines.
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


# add_region <- function(X, group, region = "who"){
#   if(!is.ic_data(X) && missing(group)){
#     stop("Please supply a grouping variable or 'ic' data.")
#   }
#
#   if(missing(groupVar)){
#     group <- attr(X, "group")
#   } else if(!is.character(group)){
#     stop("Please supply grouping variable as a character name.")
#   }
#
#   if(length(group) > 1){ group <- group[1L] }
#   if(!group %in% names(X)){ stop("Please supply a valid column name.")}
#
#   X[["region"]] <- lookup_region(X[[group]], region)
# }


# lookup_region <- function(code, region){
#   return(region_tbl[region_tbl$code == code,][[region]])
# }
