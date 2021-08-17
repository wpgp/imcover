
#' Create an ic data object
#' @description Create an ic data object from survey data on immunisation
#'   coverage.
#' @param X Any object with immunisation coverage survey data to convert into an
#'   ic object.
#' @param ... Additional parameters for names of core data elements to pass on
#'   to \code{ic_data}. See details and \code{ic_data}.
#' @param survey,sample,evidence,validity Character of the name within \code{X}
#'   which defines the core values for ic survey data.
#' @param reduce Should country-vaccine combinations which have insufficient
#'   coverage data be removed? Default is \code{TRUE}.
#' @param minSample Numeric. Minimum reported sample size needed to keep
#'   coverage records. Default is 300.
#' @return An object of class \code{ic.df} which extends \code{data.frame}-like
#'   objects with attributes to location and preserve core data elements for
#'   immunisation coverage.
#' @details \code{ic_survey} is is an extension to \code{ic_data} and is used to
#'   create a properly formed and processed dataset for immunisation coverage
#'   modelling based on survey datasets.
#'
#' The parameters \code{min} and \code{max} are only needed when \code{expand =
#' TRUE}. If \code{max} is missing, the maximum time value observed in the
#' dataset will be used.
#' @seealso \code{\link[imcover]{ic_data}}
#' @name ic_survey
#' @export
ic_survey <- function(X, ...,
                      survey = "surveyNameEnglish", sample = "Sample_Size",
                      evidence = "evidence", validity = "validity",
                      reduce = TRUE, minSample = 300){
  UseMethod('ic_survey')
}


#' @name ic_survey
#' @export
ic_survey.ic.df <- function(X, ..., dropCols = FALSE,
                            survey = "surveyNameEnglish", sample = "Sample_Size",
                            evidence = "evidence", validity = "validity",
                            reduce = TRUE, minSample = 300){
  X
}


#' @name ic_survey
#' @export
ic_survey.ic.df <- function(X, ..., dropCols = FALSE,
                            survey = "surveyNameEnglish", sample = "Sample_Size",
                            evidence = "evidence", validity = "validity",
                            reduce = TRUE, minSample = 300){

  X <- make_ic_svy(X, dropCols, survey, sample, evidence, validity,
                   reduce, minSample)

  return(X)
}


#' @name ic_survey
#' @export
ic_survey.data.frame <- function(X, ..., dropCols = FALSE,
                                 survey = "surveyNameEnglish",
                                 sample = "Sample_Size",
                                 evidence = "evidence", validity = "validity",
                                 reduce = TRUE, minSample = 300){
  if(missing(X)){
    stop("Please supply a valid dataset.")
  } else{
    X <- ic_data(X, ..., dropCols = FALSE, expand = FALSE, validate = FALSE)
  }
  stopifnot(is.ic_data(X))

  # call internal maker function
  X <- make_ic_svy(X, dropCols, survey, sample, evidence, validity,
                   reduce, minSample)

  return(X)
}


make_ic_svy <- function(X, dropCols, survey, sample, evidence, validity,
                        reduce, minSample){

  # check/confirm additional cols present
  varnames <- c(survey, evidence, validity, sample)
  chk <- !varnames %in% names(X)
  if(sum(chk) > 0){
    stop(paste(varnames[chk], collapse = " "), " not found in dataset.")
  } else{
    attr(X, "survey") <- survey
    attr(X, "sample") <- sample
    attr(X, "evidence") <- evidence
    attr(X, "validity") <- validity
  }

  # potentially missing/malformed
  X[[evidence]] <- trimws(tolower(as.character(X[[evidence]])))
  X[[validity]] <- trimws(tolower(as.character(X[[validity]])))
  # sort
  attrs <- get_attr(X, attrs = ic_core(survey = TRUE), unlist = FALSE)
  corenames <- unlist(attrs)
  X <- X[, c(corenames, setdiff(names(X), corenames))]

  # process to select preferred vacc record
  if(reduce){
    X <- survey_reduce(X, minSample)
  }

  if(dropCols){
    X <- X[, names(X) %in% get_attr(X, ic_core(survey = TRUE))]
  }

  return(X)
}


survey_adjust <- function(X, adjVacc = c("DTP", "PCV")){
  if(!is.ic_data(X)){ stop("Please supply a valid 'ic' dataset.") }
  attrs <- get_attr(X, attrs = ic_core(survey = TRUE), unlist = FALSE)
  corenames <- unlist(attrs)

  if(any(!names(attrs) %in% ic_core(survey = TRUE))){
    stop("Missing required 'ic' core attributes.")
  }

  # check vaccines (doses 1 and 3) are in the data
  vdoses <- paste0(adjVacc, rep(c(1, 3), each = length(adjVacc)))
  chk <- vdoses %in% list_vaccines(X)
  if(any(!chk)){
    stop("Missing vaccines: ", paste(vdoses[!chk], collapse = " "))
  }

  # check for missing validity info? some empty strings. drop?
  if(any(!X[[corenames['validity']]] %in% c("crude", "valid"))){
    stop(paste0("Validity must be 'crude' or 'valid'. Found: ",
                paste(unique(X[[corenames['validity']]]), collapse = ", ")),
         call. = FALSE)
  }
  # storage
  # X[["adj_factor"]] <- X[["coverage_adj"]] <- NA

  for(v in adjVacc){
    # subset vaccine to process
    vd3 <- as.data.frame(X[X[[attrs$vaccine]] == paste0(v, 3), ])
    vd1 <- as.data.frame(X[X[[attrs$vaccine]] == paste0(v, 1), ])
    vd3$adj_factor <- vd3$coverage_adj <- NULL
    vd1$adj_factor <- vd1$coverage_adj <- NULL

    # alternative to reshaping
    # vd1 <- as.data.frame(X[X[[attrs$vaccine]] == paste0(v, 1), ])
    vd1.c <- vd1[vd1[[corenames["evidence"]]] %in% c("card"), ]
    vd1.coh <- vd1[vd1[[corenames["evidence"]]] %in% c("card or history"), ]

    # find preferred data
    vd1s <- split(vd1.c, vd1.c[, corenames[c("country", "time")]], drop = TRUE)
    vd1s <- lapply(vd1s, FUN = function(vd){
      if(nrow(vd) > 1L){
        if(any(vd[[corenames["validity"]]] == "valid")){
          i <- which(vd[[corenames["validity"]]] == "valid")
          if(length(i) == 1) { vd <- vd[i, ] }
          else{ vd <- vd[which.max(vd[[corenames["sample"]]]), ] }
        } else{ # get largest sample
          vd <- vd[which.max(vd[[corenames["sample"]]]), ]
        }
      }
      return(vd)
    })
    vd1.c <- do.call(rbind, vd1s)

    vd1s <- split(vd1.coh, vd1.coh[, corenames[c("country", "time")]], drop = TRUE)
    vd1s <- lapply(vd1s, FUN = function(vd){
      if(nrow(vd) > 1L){
        if(any(vd[[corenames["validity"]]] == "valid")){
          v <- which(vd[[corenames["validity"]]] == "valid")
          if(length(v) == 1) { vd <- vd[v, ] }
          else{ vd <- vd[which.max(vd[[corenames["sample"]]]), ] }
        } else{ # get largest sample
          vd <- vd[which.max(vd[[corenames["sample"]]]), ]
        }
      }
      return(vd)
    })
    vd1.coh <- do.call(rbind, vd1s)

    # calculate adjustment factors
    vd1.coh <- merge(vd1.coh, vd1.c[, corenames[c("country","time","coverage")]],
                     by = corenames[c("country","time")],
                     suffixes = c("", ".c"))
    vd1.coh$adj_factor <- vd1.coh[[corenames["coverage"]]] / vd1.coh[[paste0(corenames["coverage"], ".c")]]


    # repeat for dose 3
    vd3.c <- vd3[vd3[[corenames["evidence"]]] %in% c("card"), ]
    vd3.coh <- vd3[vd3[[corenames["evidence"]]] %in% c("card or history"), ]

    # find preferred data
    vd3s <- split(vd3.c, vd3.c[, corenames[c("country", "time")]], drop = TRUE)
    vd3s <- lapply(vd3s, FUN = function(vd){
      if(nrow(vd) > 1L){
        if(any(vd[[corenames["validity"]]] == "valid")){
          i <- which(vd[[corenames["validity"]]] == "valid")
          if(length(i) == 1) { vd <- vd[i, ] }
          else{ vd <- vd[which.max(vd[[corenames["sample"]]]), ] }
        } else{ # get largest sample
          vd <- vd[which.max(vd[[corenames["sample"]]]), ]
        }
      }
      return(vd)
    })
    vd3.c <- do.call(rbind, vd3s)

    # calculate adjusted coverage
    vd3.c <- merge(vd3.c,
                   vd1.coh[, c(corenames[c("country","time")], "adj_factor")],
                   by = corenames[c("country", "time")])

    vd3.c$coverage_adj <- vd3.c[[corenames["coverage"]]] * vd3.c[["adj_factor"]]

    # update vd3.coh
    vd3s <- split(vd3.coh, vd3.coh[, corenames[c("country", "time")]], drop = TRUE)
    vd3s <- lapply(vd3s, FUN = function(vd){
      if(nrow(vd) > 1L){
        if(any(vd[[corenames["validity"]]] == "valid")){
          i <- which(vd[[corenames["validity"]]] == "valid")
          if(length(i) == 1) { vd <- vd[i, ] }
          else{ vd <- vd[which.max(vd[[corenames["sample"]]]), ] }
        } else{ # get largest sample
          vd <- vd[which.max(vd[[corenames["sample"]]]), ]
        }
      }
      return(vd)
    })
    vd3.coh <- do.call(rbind, vd3s)

    vd3.coh <- merge(vd3.coh,
                     vd3.c[, c(corenames[c("country","time")],
                               "adj_factor", "coverage_adj")],
                     by = corenames[c("country", "time")],
                     all.x = TRUE)

    attributes(vd3.coh)[names(attrs)] <- attrs
    class(vd3.coh) <- c("ic.df", class(vd3.coh))


    # update main dataset
    X <- X[!(X[[attrs$vaccine]] == paste0(v, 3) &
             X[[corenames["evidence"]]] %in% c("card or history")), ]
    X <- rbind(X, vd3.coh, fill = TRUE)
    X[[corenames["coverage"]]] <- ifelse(!is.na(X[["coverage_adj"]]),
                                         X[["coverage_adj"]],
                                         X[[corenames["coverage"]]])
    X <- X[do.call(order, X[, corenames[c("country", "time", "vaccine")], drop = TRUE]), ]
  }
  return(X)
}


survey_reduce <- function(X,
                          minSample = 300,
                          priority = c("card or history", "card")){
  if(!is.ic_data(X)){ stop("Please supply a valid 'ic' dataset.") }
  attrs <- get_attr(X, attrs = ic_core(survey = TRUE), unlist = FALSE)
  corenames <- unlist(attrs)

  # drop non-card or non-coh records??
  X <- as.data.frame(X)
  # drop other records
  # X <- X[!is.na(X[[attrs$sample]]) | !is.na(X[[attrs$validity]]), ]
  X <- X[(!is.na(X[[attrs$sample]]) & X[[attrs$sample]] >= minSample) |
           (!is.na(X[[attrs$validity]]) & X[[attrs$validity]] == "valid"), ]

  # create country x time x vaccine sets to process
  xs <- split(X, X[, corenames[c("country", "time", "vaccine")]], drop = TRUE)
  xs <- lapply(xs, FUN = function(x){
    # x <- xs[[i]]
    if(nrow(x) == 1L){
      return(x)
    } else{ # multiple rows
      # find preferred records
      valid_f <- factor(x[[attrs$validity]], levels = c("valid", "crude"))
      evid_f <- factor(x[[attrs$evidence]], levels = priority)
      x <- x[order(valid_f, evid_f, x[[attrs$sample]], decreasing = c(FALSE, FALSE, TRUE)), ]
      return(x[1L, ])
    } # end else
  })
  X <- do.call(rbind.data.frame, xs)
  row.names(X) <- seq(nrow(X))

  attributes(X)[names(attrs)] <- attrs
  # attr(X, "class") <- c("ic.svy", "ic.df", class(X))
  attr(X, "class") <- c("ic.df", class(X))
  stopifnot(is.ic_data(X))

  return(X)
}


#' #' Find survey groups
#' mark_survey <- function(x){
#'   if(!is.ic_data(x)){ stop("Please supply a valid 'ic' dataset.") }
#'   vars <- get_attr(x, c("country", "survey", "time", "vaccine", "coverage", "sample"))
#'
#'   x <- as.data.frame(x)
#'   x$N <- ave(x[[vars["coverage"]]], x[, vars[c("country", "survey", "time")]], FUN = length)
#'
#'   x <- x[do.call(order, x[, vars[c("country", "survey", "time")]]), ]
#'
#'   x <- within(x, n <- ave(get_attr(x, "coverage"),
#'                           as.list(c(vars, get_attr(x, "evidence"))),
#'                           FUN = seq_along))
#'   return(x)
#' }

