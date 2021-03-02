

#' Create an ic data object
#' @description Create an ic data object from survey data on immunisation
#'   coverage.
#' @param X Any object with immunisation coverage survey data to convert into an
#'   ic object.
#' @param ... Additional parameters for names of core data elements to pass on
#'   to \code{ic_data}. See details and \code{ic_data}.
#' @param survey,sample,evidence,validity Character of the name within \code{X}
#'   which defines the core value for ic survey data.
#' @param reduce Should group-vaccine combinations which have no coverage data
#'   be removed? Default is \code{TRUE}.
#' @param minSample Numeric. Minimum reported sample size. Default is 300.
#' @param expand Should additional empty rows be added to expand the ic data to
#'   all possible group x vaccine x time observation? Default is \code{FALSE}.
#' @param min Start of the potential time period. Default is 1999.
#' @param max End of the potential time period.
#' @return An object of class \code{ic.svy} which extends \code{ic.df} and
#'   \code{data.frame}-like objects with attributes to location and preserve
#'   core data elements for immunisation coverage.
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
ic_survey <- function(X, ..., dropCols = FALSE,
                      survey = "survey_name", sample = "sample_size",
                      evidence = "evidence", validity = "validity",
                      reduce = TRUE, minSample = 300,
                      expand = TRUE, min = 1999, max, na.remove = TRUE,
                      biasAdjust = TRUE, adjVacc = c("DTP", "PCV")){
  UseMethod('ic_survey')
}


#' @name ic_survey
#' @export
ic_survey.ic.svy <- function(X, ..., dropCols = FALSE,
                             survey = "survey_name", sample = "sample_size",
                             evidence = "evidence", validity = "validity",
                             reduce = TRUE, minSample = 300,
                             expand = TRUE, min = 1999, max, na.remove = TRUE,
                             biasAdjust = TRUE, adjVacc = c("DTP", "PCV")){
  X
}


#' @name ic_survey
#' @export
ic_survey.ic.df <- function(X, ..., dropCols = FALSE,
                            survey = "survey_name", sample = "sample_size",
                            evidence = "evidence", validity = "validity",
                            reduce = TRUE, minSample = 300,
                            expand = TRUE, min = 1999, max, na.remove = TRUE,
                            biasAdjust = TRUE, adjVacc = c("DTP", "PCV")){
  X # add processing steps for bias and sample corrections
}


#' @name ic_survey
#' @export
ic_survey.data.frame <- function(X, ..., dropCols = FALSE,
                                 survey = "survey_name", sample = "sample_size",
                                 evidence = "evidence", validity = "validity",
                                 reduce = TRUE, minSample = 300,
                                 expand = TRUE, min = 1999, max, na.remove = TRUE,
                                 biasAdjust = TRUE, adjVacc = c("DTP", "PCV")){
  if(missing(X)){
    stop("Please supply a valid dataset.")
  } else{
    X <- ic_data(X, ..., dropCols = FALSE, expand = FALSE, validate = FALSE)
  }
  stopifnot(is.ic_data(X))

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

  # class(X) <- c("ic.svy", class(X))
  # potentially missing/malformed
  X[[evidence]] <- trimws(tolower(as.character(X[[evidence]])))
  X[[validity]] <- trimws(tolower(as.character(X[[validity]])))
  # sort
  attrs <- get_attr(X, attrs = ic_core(survey = TRUE), unlist = FALSE)
  corenames <- unlist(attrs)
  X <- X[, c(corenames, setdiff(names(X), corenames))]

  # correct dose 3 recall
  if(biasAdjust){
    X <- survey_adjust(X, adjVacc)
  }

  # process to select preferred vacc record
  if(reduce){
    X <- survey_reduce(X, minSample)
  }

  # add in empty rows
  if(expand){
    if(missing(max)) max <- base::max(X[[attrs$time]], na.rm = TRUE)
    X <- ic_expand(X, min = min, max = max, na.remove = na.remove)
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
    vd1s <- split(vd1.c, vd1.c[, corenames[c("group", "time")]], drop = TRUE)
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

    vd1s <- split(vd1.coh, vd1.coh[, corenames[c("group", "time")]], drop = TRUE)
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
    vd1.coh <- merge(vd1.coh, vd1.c[, corenames[c("group","time","coverage")]],
                     by = corenames[c("group","time")],
                     suffixes = c("", ".c"))
    vd1.coh$adj_factor <- vd1.coh[[corenames["coverage"]]] / vd1.coh[[paste0(corenames["coverage"], ".c")]]


    # repeat for dose 3
    vd3.c <- vd3[vd3[[corenames["evidence"]]] %in% c("card"), ]
    vd3.coh <- vd3[vd3[[corenames["evidence"]]] %in% c("card or history"), ]

    # find preferred data
    vd3s <- split(vd3.c, vd3.c[, corenames[c("group", "time")]], drop = TRUE)
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
                   vd1.coh[, c(corenames[c("group","time")], "adj_factor")],
                   by = corenames[c("group", "time")])

    vd3.c$coverage_adj <- vd3.c[[corenames["coverage"]]] * vd3.c[["adj_factor"]]

    # update vd3.coh
    vd3s <- split(vd3.coh, vd3.coh[, corenames[c("group", "time")]], drop = TRUE)
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
                     vd3.c[, c(corenames[c("group","time")],
                               "adj_factor", "coverage_adj")],
                     by = corenames[c("group", "time")],
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
    X <- X[do.call(order, X[, corenames[c("group", "time", "vaccine")], drop = TRUE]), ]

    # # reshape to wide - collapsing validity information
    # vd1 <- reshape(vd1, direction = "wide",
    #                idvar = corenames[c("group","survey","time",
    #                                    "vaccine","evidence")],
    #                timevar = corenames["validity"])
    # # find preferred dose 1 data
    # miss <- setdiff(paste(attrs$coverage, c("valid", "crude"), sep = "."),
    #                 names(vd1))
    # if(length(miss) > 1){
    #   stop("Missing validity.")
    # } else if(length(miss) == 1){
    #   vd1[[miss]] <- NA
    # }
    # prefer <- paste(attrs$coverage, "valid", sep = ".")
    #
    # vd1[[attrs$validity]] <- ifelse(is.na(vd1[[prefer]]),
    #                                 "crude", "valid")
    #
    # vd1[[attrs$coverage]] <- ifelse(is.na(vd1[[prefer]]),
    #                                 vd1[[paste(attrs$coverage, "crude", sep = ".")]],
    #                                 vd1[[paste(attrs$coverage, "valid", sep = ".")]])
    #
    # vd1[[attrs$sample]] <- ifelse(is.na(vd1[[prefer]]),
    #                               vd1[[paste(attrs$sample, "crude", sep = ".")]],
    #                               vd1[[paste(attrs$sample, "valid", sep = ".")]])
    #
    # # merge c and coh
    # vd1.coh <- merge(vd1[vd1[[attrs$evidence]] == "card or history", corenames],
    #                  vd1[vd1[[attrs$evidence]] == "card", corenames],
    #                  by = corenames[c('group','time','survey')],
    #                  all.x = TRUE, suffixes = c(".coh", ".c"))
    # vd1 <- vd1[!vd1[[attrs$evidence]] %in% c("card", "card or history"), corenames]
    # # calculate adjustment factor
    # vd1.coh$adj_factor <- vd1.coh[, paste0(attrs$coverage, ".coh")] / vd1.coh[, paste0(attrs$coverage, ".c")]
    #
    # # process dose 3
    # # reshape to wide - collapsing validity information
    # vd3 <- reshape(vd3, direction = "wide",
    #                idvar = corenames[c("group","survey","time",
    #                                    "vaccine","evidence")],
    #                timevar = corenames["validity"])
    # # find preferred dose 3 data
    # miss <- setdiff(paste(attrs$coverage, c("valid", "crude"), sep = "."),
    #                 names(vd3))
    # if(length(miss) > 1){
    #   stop("Missing validity.")
    # } else if(length(miss) == 1){
    #   vd3[[miss]] <- NA
    # }
    # prefer <- paste(attrs$coverage, "valid", sep = ".")
    #
    # vd3[[attrs$validity]] <- ifelse(is.na(vd3[[prefer]]),
    #                                 "crude", "valid")
    #
    # vd3[[attrs$coverage]] <- ifelse(is.na(vd3[[prefer]]),
    #                                 vd3[[paste(attrs$coverage, "crude", sep = ".")]],
    #                                 vd3[[paste(attrs$coverage, "valid", sep = ".")]])
    #
    # vd3[[attrs$sample]] <- ifelse(is.na(vd3[[prefer]]),
    #                               vd3[[paste(attrs$sample, "crude", sep = ".")]],
    #                               vd3[[paste(attrs$sample, "valid", sep = ".")]])
    #
    # # adjust dose 3 'card' only data
    # vd_adj <- merge(vd3[vd3[[attrs$evidence]] == "card", corenames],
    #                 vd1.coh[, c(corenames[c('group','time','survey')], "adj_factor")],
    #                 by = corenames[c('group','time','survey')],
    #                 all.x = TRUE, suffixes = c(".vd3", ".vd1"))
    # vd_adj$coverage_adj <- vd_adj[[attrs$coverage]] * vd_adj$adj_factor
#
#     # update dose 3 'card or history' data
#     X <- merge(X,
#                vd_adj[, c(corenames[c('group','time','vaccine',
#                                       'survey','evidence','validity')],
#                           "adj_factor", "coverage_adj")],
#                by = corenames[c('group','time','vaccine','survey','evidence','validity')],
#                all.x = TRUE, suffixes = c("", ".new"))
#
#     X[["coverage_adj"]] <- ifelse(is.na(X[["coverage_adj"]]),
#                                   X[["coverage_adj.new"]],
#                                   X[["coverage_adj"]])
#
#     X[["adj_factor"]] <- ifelse(is.na(X[["adj_factor"]]),
#                                 X[["adj_factor.new"]],
#                                 X[["adj_factor"]])
#     # drop temporary columns
#     X[["coverage_adj.new"]] <- X[["adj_factor.new"]] <- NULL
  }
  return(X)
}


survey_reduce <- function(X, minSample = 300){
  if(!is.ic_data(X)){ stop("Please supply a valid 'ic' dataset.") }
  attrs <- get_attr(X, attrs = ic_core(survey = TRUE), unlist = FALSE)
  corenames <- unlist(attrs)

  # drop non-card or non-coh records??
  X <- as.data.frame(X)
  # store processing code (aka 'samseen')
  # X[["pcode"]] <- NA

  # create group x time x vaccine sets to process
  xs <- split(X, X[, corenames[c("group", "time", "vaccine")]], drop = TRUE)
  xs <- lapply(xs, FUN = function(x){
    # x <- xs[[i]]
    if(nrow(x) == 1L){
      if((!is.na(x[[attrs$sample]]) && x[[attrs$sample]] > minSample) || x[[attrs$validity]] == "valid"){
        # x[["pcode"]] <-
        return(x)
      } else{
        return(NULL) # skipped in binding
      }
    } else{
      # find preferred records
      # xc <- x[grepl("c or h|card or history", x[[attrs$evidence]]), ]
      xc <- x[c(grep("card or history", x[[attrs$evidence]], value = FALSE),
                grep("c or h", x[[attrs$evidence]], value = FALSE)), ]
      if(nrow(xc) >= 1L){
        xc <- xc[order(xc[[attrs$sample]], decreasing = TRUE), ]
        return(xc[1L, ])
      } else{
        # xc <- x[grepl("c|card", x[[attrs$evidence]]), ]
        xc <- x[c(grep("card", x[[attrs$evidence]], value = FALSE),
                  grep("c", x[[attrs$evidence]], value = FALSE)), ]
        xc <- xc[(!is.na(xc[[attrs$sample]]) && xc[[attrs$sample]] > minSample) || xc[[attrs$validity]] == "valid", ]
        if(nrow(xc) >= 1L){
          xc <- xc[order(xc[[attrs$sample]], decreasing = TRUE), ]
          return(xc[1L, ])
        } else{
          return(NULL)
        }
      }
    }
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
#'   vars <- get_attr(x, c("group", "survey", "time", "vaccine", "coverage", "sample"))
#'
#'   x <- as.data.frame(x)
#'   x$N <- ave(x[[vars["coverage"]]], x[, vars[c("group", "survey", "time")]], FUN = length)
#'
#'   x <- x[do.call(order, x[, vars[c("group", "survey", "time")]]), ]
#'
#'   x <- within(x, n <- ave(get_attr(x, "coverage"),
#'                           as.list(c(vars, get_attr(x, "evidence"))),
#'                           FUN = seq_along))
#'   return(x)
#' }

