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
#' @name ic_core
#' @export
ic_core <- function(survey = FALSE){
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
