

#' Create an ic data object
#'
#' Create an ic data object, which extends \code{data.frame}-like objects by
#' defining consistent attributes for required data elements.
#' @param X An \code{R} object to process and convert to ic data.
#' @param group,time,vaccine,coverage,dose,population Character of the name
#'   within \code{X} which defines the core value for ic data.
#' @param dropCols Should other columns in \code{X} be dropped if they are not
#'   core attributes? Default is \code{FALSE} to retain all data from \code{X}.
#' @param expand Should additional empty rows be added to expand the ic data to
#'   all possible group x vaccine x time observation? Default is \code{FALSE}.
#' @param validate Should basic ic data checking be performed? Default is
#'   \code{FALSE}.
#' @param ... Additional arguments to pass to \code{ic_expand}.
#' @return An object of class \code{ic.df} which extends \code{data.frame}-like
#'   objects with attributes to location and preserve core data elements for
#'   immunisation coverage.
#' @details \code{ic_data} is the core function to create a properly formed and
#'   processed dataset for immunisation coverage modelling. In particular it
#'   requires some core data elements are present:
#' \itemize{
#'   \item{'group'}{ Aggregate grouping for data records (e.g. country). Group
#'   can be defined by more than column for a nested hierarchy of units (largest
#'   to finest grouping). Default name is 'iso3countrycode'.}
#'   \item{'time'}{ Defines the time period of immunisation records, typically an
#'   integer year. Default name is 'year'.}
#'   \item{'vaccine'}{ Code to identify the vaccine records (e.g. 'DTP1').
#'   Default name is 'vaccine_name'.}
#'   \item{'coverage'}{ Pre-calculated coverage percentage. Default name is
#'   'percentcoverage'.}
#' }
#'   If 'coverage' is not included, then two other elements are required to be
#'   specified in order for percent coverage to be calculated. Else, these are
#'   optional elements for ic data.
#' \itemize{
#'   \item{'dose'}{ Number of vaccine doses administered. Default name is
#'   'dosesadministered'.}
#'   \item{'population'}{ Total target population for the vaccine. Default name
#'   is 'targetgroup'.}
#' }
#' @seealso \code{ic_expand}, \code{ic_validate}
#' @name ic_data
#' @export
ic_data <- function(X, group = 'iso3countrycode', time = 'year',
                    vaccine = 'vaccine_name', coverage = 'percentcoverage',
                    dose = 'dosesadministered', population = 'targetgroup',
                    dropCols = FALSE, expand = FALSE, validate = FALSE, ...){
  UseMethod("ic_data")
}


#' @name ic_data
#' @export
ic_data.ic.df <- function(X, group = 'iso3countrycode', time = 'year',
                          vaccine = 'vaccine_name', coverage = 'percentcoverage',
                          dose = 'dosesadministered', population = 'targetgroup',
                          dropCols = FALSE, expand = FALSE, validate = FALSE, ...){
  X
}


#' @name ic_data
#' @export
ic_data.data.frame <- function(X,
                               group = 'iso3countrycode', time = 'year',
                               vaccine = 'vaccine_name', coverage = 'percentcoverage',
                               dose = 'dosesadministered', population = 'targetgroup',
                               dropCols = FALSE, expand = FALSE, validate = FALSE, ...){

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

  varnames <- c(group, time, vaccine)
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
  attr(X, "group") <- group
  attr(X, "time") <- time
  attr(X, "vaccine") <- vaccine
  attr(X, "coverage") <- coverage
  attr(X, "dose") <- dose
  attr(X, "population") <- population

  if(expand){
    X <- ic_expand(X, ...)
  }

  # if(validate){
  #   X <- ic_validate(X)
  # }

  return(X)
}


#' Convert foreign object to an ic data object
#'
#' @description Functions to check if an object is an ic data object, or coerce
#'   it to one if possible.
#' @param object Any \code{R} object.
#' @param ... Elements to be combined into an \code{ic.df} object.
#' @details Core \code{ic.df} data elements are expected. Specifically, data for
#'   'group', 'time', 'vaccine', (optionally 'dose', 'population'), and
#'   'coverage' are required in this order. The names of these elements from
#'   \code{...} are passed to \code{ic_data}.
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

  if(length(varnames) < 4){
    stop("Not enough columns. Expecting data for: group, time, vaccine,
         (doses, target population,) and coverage.")
  }

  # set attributes
  class(x) <- list("ic.df", class(x))
  if(length(varnames) == 4){
    attributes(x)[c("group", "time", "vaccine", "coverage")] <- varnames[1:4]
  } else{
    attributes(x)[c("group", "time", "vaccine",
                    "dose", "population", "coverage")] <- varnames[1:6]
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


#' Expand ic data to include all potential years.
#'
#' Expand the time dimension of an ic dataset to include the full set of
#' potential observations.
#' @param X Object of class \code{ic.df} to expand.
#' @param min Start of the potential time period. Default is 1999.
#' @param max End of the potential time period.
#' @param na.remove Should group-vaccine combinations which have no coverage
#'   data be removed? Default is \code{TRUE}
#' @details If \code{max} is missing, the maximum time value observed in the
#'   dataset will be used.
#' @name ic_expand
#' @export
ic_expand <- function(X, min = 1999, max, na.remove = TRUE){
  if(!is.ic_data(X)){
    stop("Please provide valid 'ic' data.")
  }
  # check years
  if(missing(max) || is.na(max)){
    max <- base::max(X[[attr(X, "time")]], na.rm = TRUE)
  }

  if(min > max){
    stop("Invalid timespan.")
  }
  years <- min:max

  # expand groups
  df <- unique(X[, c(attr(X, "group"), attr(X, "vaccine")), drop = TRUE])
  times <- rep(years, times=nrow(df))
  df <- df[rep(seq_len(nrow(df)), each=length(years)), ]
  df[[attr(X, "time")]] <- times
  # add NAs
  df <- merge(X, df, # merge.ic.df
              by = c(attr(X, "group"), attr(X, "time"), attr(X, "vaccine")),
              all.y = TRUE, sort = FALSE, attr.x = TRUE)

  # drop group x vacc where all years = NA
  if(na.remove){
    splits <- split(df,
                    f = df[, c(attr(X, "group"), attr(X, "vaccine")), drop = TRUE],
                    drop = TRUE)
    empty <- lapply(splits,
                    FUN = function(i){ all(is.na(i[[attr(X, "coverage")]])) })
    empty <- unlist(empty, use.names = FALSE)

    df <- do.call(rbind.data.frame, splits[which(!empty)])
  }

  # sort
  sort_list <- c(attr(X, "group"), attr(X, "time"), attr(X, "vaccine"))
  df <- df[do.call(order, df[ , match(sort_list, names(df))]),]
  row.names(df) <- seq(nrow(df))

  # df <- do.call("ic_data", c(list(X=df), attributes(X)))
  return(df)
}


# ic_update <- function(X, compare, validate = TRUE){
#   if(!is.ic_data(X)){
#     stop("Please provide valid 'ic' data.")
#   }
#
# }


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


# get_yovi <- function(X){
#   # return full yovi table
# }


#' @name ic_update_pop
#' @export
ic_update_pop <- function(X, pop = "wpp", validate = TRUE){
  if(!is.ic_data(X)){
    stop("Please provide valid 'ic' data.")
  }

  TRUE
}


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



