
#' Create an ic data object
#'
#' Create an ic data object, which extends \code{data.frame}-like objects by
#' defining consistent attributes for required data elements.
#' @param X An \code{R} object to process and convert to ic data.
#' @param region,country,time,vaccine,coverage,dose,population Character of the
#'   name within \code{X} which defines the core value for ic data.
#' @param dropCols Should other columns in \code{X} be dropped if they are not
#'   core attributes? Default is \code{FALSE} to retain all data from \code{X}.
#' @param ... Additional arguments. Not currently used.
#' @return An object of class \code{ic.df} which extends \code{data.frame}-like
#'   objects with attributes to location and preserve core data elements for
#'   immunisation coverage.
#' @details \code{ic_data} is the core function to create a properly formed and
#'   processed dataset for immunisation coverage modelling. In particular it
#'   requires some core data elements are present:
#' \itemize{
#'   \item{'region'}{ Larger grouping. Used to apply separate models.}
#'   \item{'country'}{ Aggregate grouping for data records (e.g. country). Country
#'   can be defined by more than column for a nested hierarchy of units (largest
#'   to finest grouping). Default name is 'code' referring to an ISO3 code.}
#'   \item{'time'}{ Defines the time period of immunisation records, typically
#'   an integer year. Default name is 'year'.}
#'   \item{'vaccine'}{ Code to identify the vaccine records (e.g. 'DTP1').
#'   Default name is 'antigen'.}
#'   \item{'coverage'}{ Pre-calculated coverage percentage. Default name is
#'   'coverage'.}
#'   \item{'source'}{ Character identifying the source of the data. Default is
#'   'coverage_category'.}
#' }
#'   If 'coverage' is not included, then two other elements are required to be
#'   specified in order for percent coverage to be calculated. Else, these are
#'   optional elements for ic data.
#' \itemize{
#'   \item{'dose'}{ Number of vaccine doses administered. Default name is
#'   'doses'.}
#'   \item{'population'}{ Total target population for the vaccine. Default name
#'   is 'target_number.}
#' }
#'
#' @return An object of type \code{ic.df}.
#'
#' @seealso \code{\link[imcover]{ic_expand}}, \code{ic_validate}
#' @name ic_data
#' @export
ic_data <- function(X, region = 'region', country = 'code', time = 'year',
                    vaccine = 'antigen', coverage = 'coverage',
                    source = 'coverage_category', dose = 'doses',
                    population = 'target_number', dropCols = FALSE, ...){
  UseMethod("ic_data")
}


#' @name ic_data
#' @export
ic_data.ic.df <- function(X, region = 'region', country = 'code', time = 'year',
                          vaccine = 'antigen', coverage = 'coverage',
                          source = 'coverage_category', dose = 'doses',
                          population = 'target_number', dropCols = FALSE, ...){
  X
}


#' @name ic_data
#' @export
ic_data.data.frame <- function(X, region = 'region', country = 'code', time = 'year',
                               vaccine = 'antigen', coverage = 'coverage',
                               source = 'coverage_category', dose = 'doses',
                               population = 'target_number', dropCols = FALSE, ...){

  # check data types
  if(missing(X)){
    stop("Please supply a valid dataset.")
  }

  if(inherits(X, "list")){
    if(is.null(names(X))){ stop("Please supply valid data names.") }
    X <- data.frame(X, stringsAsFactors = FALSE)
  }

  # check names
  if(any(lengths(list(time, vaccine, coverage, dose, population)) > 1)){
    stop("Please provide valid column names.")
  }

  varnames <- c(region, country, time, vaccine, source)
  chk <- !varnames %in% names(X)
  if(sum(chk) > 0){
    stop(paste(varnames[chk], collapse = " "), " not found in dataset.")
  }
  # remove factors
  for(i in varnames){
    if(is.factor(X[[i]])) X[[i]] <- as.character(X[[i]])
  }
  # clean vax names
  X[[vaccine]] <- trimws(toupper(as.character(X[[vaccine]])))

  # get vaccine coverage
  if(!coverage %in% names(X)){
    if(dose %in% names(X) && population %in% names(X)){
      X[[coverage]] <- X[[dose]] / X[[population]] * 100
    } else{
      stop("Provide coverage data or doses and target population.")
    }
  }
  if(!dose %in% names(X)){ dose <- NULL }
  if(!population %in% names(X)){ population <- NULL }

  keepnames <- c(varnames, dose, population, coverage)
  if(dropCols){
    X <- X[, keepnames]
  } else{
    othernames <- names(X)[!names(X) %in% keepnames]
    X <- X[, c(keepnames, othernames)]
  }

  # set attributes
  class(X) <- list("ic.df", class(X))
  attr(X, "region") <- region
  attr(X, "country") <- country
  attr(X, "time") <- time
  attr(X, "vaccine") <- vaccine
  attr(X, "coverage") <- coverage
  attr(X, "source") <- source
  attr(X, "dose") <- dose
  attr(X, "population") <- population

  return(X)
}


#' Convert foreign object to an ic data object
#'
#' @description Functions to check if an object is an ic data object, or coerce
#'   it to one if possible.
#' @param object Any \code{R} object.
#' @param ... Elements to be combined into an \code{ic.df} object.
#' @details Core \code{ic.df} data elements are expected. Specifically, data for
#'   'region', 'country', 'time', 'vaccine', 'coverage', and 'source' (optionally
#'   'dose', 'population') are required in this order. The names of these
#'   elements from \code{...} are passed to \code{ic_data}.
#'
#' @return An object of type \code{ic.df}.
#'
#' @seealso \code{ic_data}
#' @name as.ic_data
#' @export
as.ic_data <- function(...){
  x <- list(...)

  if(length(x) == 1L && inherits(x[[1L]], c("data.frame", "tbl_df"))){
    x <- x[[1L]]
  } else{
    x <- data.frame(x, stringsAsFactors = FALSE)
  }

  varnames <- names(x)
  if(any(is.null(varnames))){ stop("Missing data names not allowed.") }

  if(length(varnames) < 6){
    stop("Not enough columns. Expecting data for: country, time, vaccine,
         (doses, target population,) and coverage.")
  }

  # set attributes
  class(x) <- list("ic.df", class(x))
  if(length(varnames) == 6){
    attributes(x)[c("region", "country", "time",
                    "vaccine", "coverage", "source")] <- varnames[1:6]
  } else{
    attributes(x)[c("region", "country", "time", "vaccine",
                    "coverage", "source", "dose",
                    "population")] <- varnames[1:8]
  }

  return(x)
}


#' @name as.ic_data
#' @export
is.ic_data <- function(object){
  return(inherits(object, "ic.df"))
}


#' Access ic data
#'
#' @description Operators to extract or replace parts of ic data.
#' @param x Object of class \code{ic.df} from which to extract element(s) or in
#'   which to replace element(s).
#' @param i,j,... Indices specifying the elements to extract or replace. Numeric
#'   or character vector or empty (missing) or NULL.
#' @param drop Can core ic attribute columns be modified? Default is
#'   \code{FALSE}.
#' @details \code{[.ic.df} will generally follow the behaviour of
#'   \code{data.frame}; however, when \code{drop = FALSE} the core 'ic'
#'   attributes will always be preserved and will always return an object of
#'   class \code{ic.df}.
#' @seealso \code{ic_data}
#' @name ic.df
#' @export
"[.ic.df" <- function(x, i, j, ..., drop = FALSE){
  mdrop <- missing(drop)
  nargs <- nargs() - !mdrop

  attrs <- attributes(x)
  cattrs <- attrs[names(attrs) %in% ic_core()]

  coredf <- as.data.frame(x)[, unlist(cattrs)]

  if(!missing(i) && nargs > 2){ # x[3:4,]
    if(is.character(i)) i <- match(i, row.names(x))
    coredf <- coredf[i,]
  }

  class(x) <- setdiff(class(x), "ic.df") # drop to data.frame
  x <- if(missing(j)){
    if(nargs == 2) {
      x[i]
    }
    else x[i, , drop = drop]
  } else{
    x[i, j, drop = drop]
  }

  if(!drop){ # prevent loss of core ic attributes
    miss <- setdiff(unlist(cattrs), names(x))
    if(length(miss) > 0){
      x[miss] <- coredf[, miss]
    }
    x <- x[, c(unlist(cattrs), setdiff(names(x), unlist(cattrs)))]
    # remake ic data - else return unclass()
    cattrs[!cattrs %in% names(x)] <- NULL
    x <- do.call("ic_data", c(list(X=x), cattrs))

    # check/confirm additional cols present
    varnames <- ic_core(survey = TRUE)
    attributes(x)[varnames] <- attrs[varnames]
  }

  return(x)
}


#' #' @export
#' "[<-.ic.df" <- function(x, i, j, value){
#'   attrs <- get_attr(x, ic_core(), unlist = F)
#'   class(x) <- setdiff(class(x), "ic.df")
#'
#'   if(!(i %in% unlist(attrs) && is.null(value))){
#'     x <- structure(NextMethod(),
#'                    class = c("ic.df", setdiff(class(x), "ic.df")))
#'   } else{
#'     stop("Cannot drop core 'ic' data columns.", call. = FALSE, )
#'   }
#'
#'   return(x)
#' }


#' @export
"[[<-.ic.df" <- function(x, i, value){
  attrs <- get_attr(x, ic_core(), unlist = TRUE)

  if(!(i %in% attrs && is.null(value))){
    x <- structure(NextMethod(),
                   class = c("ic.df", setdiff(class(x), "ic.df")))
  } else{
    stop("Cannot drop core 'ic' data columns.", call. = FALSE, )
  }

  return(x)
}


#' @export
"[[.ic.df" <- function(x, ...){
  class(x) <- setdiff(class(x), "ic.df") # drop
  x[[...]]
}


#' @export
"$<-.ic.df" <- function(x, name, value){
  attrs <- get_attr(x, ic_core(), unlist = TRUE)

  if(name %in% attrs && is.null(value)){
    stop("Cannot drop core 'ic' data columns.", call. = FALSE, )
  } else{
    x[[name]] <- value
  }
  return(x)
}


#' @export
'names<-.ic.df' <- function(x, value){
  chknames(value)
  attrs <- get_attr(x, ic_core(survey = TRUE))
  oldc <- setdiff(class(x), "data.frame")
  oldnms <- which(names(x) %in% attrs)

  y <- as.data.frame(x)
  names(y) <- value

  attributes(y)[names(attrs)] <- names(y)[oldnms]
  class(y) <- c(oldc, class(y))

  y
}


#' @export
'names.ic.df' <- function(x){
  names(as.data.frame(x))
}


chknames <- function(x) {
  if (!identical(x, make.names(x)))
    warning("Found potentially invalid names.")
}




