

#' Create an ic data object
#' @description Create an ic data object from survey data on immunisation
#'   coverage.
#' @param X Any object with immunisation coverage survey data to convert into an
#'   ic object.
#' @param ... Additional parameters for names of core data elements to pass on
#'   to \code{ic_data}. See details and \code{ic_data}.
#' @param evidence Character within \code{X} or a vector of ...
#' @details
#' @seealso \code{ic_data}
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
  attr(X, "class") <- c("ic.svy", class(X))
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
    X <- ic_expand(X, min = min, max = max, na.remove = TRUE)
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
                paste(unique(X[[validity]]), collapse = ", ")),
         call. = FALSE)
  }
  # storage
  X[["adj_factor"]] <- X[["coverage_adj"]] <- NA

  for(v in adjVacc){
    # subset vaccine to process
    vd3 <- as.data.frame(X[X[[attrs$vaccine]] == paste0(v, 3), ])
    vd1 <- as.data.frame(X[X[[attrs$vaccine]] == paste0(v, 1), ])

    # reshape to wide - collapsing validity information
    vd1 <- reshape(vd1, direction = "wide",
                   idvar = corenames[c("group","survey","time",
                                       "vaccine","evidence")],
                   timevar = corenames["validity"])
    # find preferred dose 1 data
    miss <- setdiff(paste(attrs$coverage, c("valid", "crude"), sep = "."),
                    names(vd1))
    if(length(miss) > 1){
      stop("Missing validity.")
    } else if(length(miss) == 1){
      vd1[[miss]] <- NA
    }
    prefer <- paste(attrs$coverage, "valid", sep = ".")

    vd1[[attrs$validity]] <- ifelse(is.na(vd1[[prefer]]),
                                    "crude", "valid")

    vd1[[attrs$coverage]] <- ifelse(is.na(vd1[[prefer]]),
                                    vd1[[paste(attrs$coverage, "crude", sep = ".")]],
                                    vd1[[paste(attrs$coverage, "valid", sep = ".")]])

    vd1[[attrs$sample]] <- ifelse(is.na(vd1[[prefer]]),
                                  vd1[[paste(attrs$sample, "crude", sep = ".")]],
                                  vd1[[paste(attrs$sample, "valid", sep = ".")]])

    # merge c and coh
    vd1.coh <- merge(vd1[vd1[[attrs$evidence]] == "card or history", corenames],
                     vd1[vd1[[attrs$evidence]] == "card", corenames],
                     by = corenames[c('group','time','survey')],
                     all.x = TRUE, suffixes = c(".coh", ".c"))
    vd1 <- vd1[!vd1[[attrs$evidence]] %in% c("card", "card or history"), corenames]
    # calculate adjustment factor
    vd1.coh$adj_factor <- vd1.coh[, paste0(attrs$coverage, ".coh")] / vd1.coh[, paste0(attrs$coverage, ".c")]

    # process dose 3
    # reshape to wide - collapsing validity information
    vd3 <- reshape(vd3, direction = "wide",
                   idvar = corenames[c("group","survey","time",
                                       "vaccine","evidence")],
                   timevar = corenames["validity"])
    # find preferred dose 3 data
    miss <- setdiff(paste(attrs$coverage, c("valid", "crude"), sep = "."),
                    names(vd3))
    if(length(miss) > 1){
      stop("Missing validity.")
    } else if(length(miss) == 1){
      vd3[[miss]] <- NA
    }
    prefer <- paste(attrs$coverage, "valid", sep = ".")

    vd3[[attrs$validity]] <- ifelse(is.na(vd3[[prefer]]),
                                    "crude", "valid")

    vd3[[attrs$coverage]] <- ifelse(is.na(vd3[[prefer]]),
                                    vd3[[paste(attrs$coverage, "crude", sep = ".")]],
                                    vd3[[paste(attrs$coverage, "valid", sep = ".")]])

    vd3[[attrs$sample]] <- ifelse(is.na(vd3[[prefer]]),
                                  vd3[[paste(attrs$sample, "crude", sep = ".")]],
                                  vd3[[paste(attrs$sample, "valid", sep = ".")]])

    # adjust dose 3 'card' only data
    vd_adj <- merge(vd3[vd3[[attrs$evidence]] == "card", corenames],
                    vd1.coh[, c(corenames[c('group','time','survey')], "adj_factor")],
                    by = corenames[c('group','time','survey')],
                    all.x = TRUE, suffixes = c(".vd3", ".vd1"))
    vd_adj$coverage_adj <- vd_adj[[attrs$coverage]] * vd_adj$adj_factor

    # update dose 3 'card or history' data
    X <- merge(X,
               vd_adj[, c(corenames[c('group','time','vaccine',
                                      'survey','evidence','validity')],
                          "adj_factor", "coverage_adj")],
               by = corenames[c('group','time','vaccine','survey','evidence','validity')],
               all.x = TRUE, suffixes = c("", ".new"))

    X[["coverage_adj"]] <- ifelse(is.na(X[["coverage_adj"]]),
                                  X[["coverage_adj.new"]],
                                  X[["coverage_adj"]])

    X[["adj_factor"]] <- ifelse(is.na(X[["adj_factor"]]),
                                X[["adj_factor.new"]],
                                X[["adj_factor"]])
    # drop temporary columns
    X[["coverage_adj.new"]] <- X[["adj_factor.new"]] <- NULL
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
    x <- xs[[i]]
    if(nrow(x) == 1L){
      if((!is.na(x[[attrs$sample]]) && x[[attrs$sample]] > minSample) || x[[attrs$validity]] == "valid"){
        # x[["pcode"]] <-
        return(x)
      } else{
        return(NULL) # skipped in binding
      }
    } else{
      # find preferred records
      xc <- x[grepl("c or h|card or history", x[[attrs$evidence]]), ]
      if(nrow(xc) >= 1L){
        xc <- xc[order(xc[[attrs$sample]], decreasing = TRUE), ]
        return(xc[1L, ])
      } else{
        xc <- x[grepl("c|card", x[[attrs$evidence]]), ]
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
  attr(X, "class") <- c("ic.svy", "ic.df", class(X))
  stopifnot(is.ic_data(X))

  return(X)
}


#' #' Find survey groups
#' mark_survey <- function(x){
#'   vars <- get_attr(x, c("group", "survey", "time", "vaccine"))
#'   # x <- within(x, N <- ave(get_attr(x, "coverage"), as.list(vars), FUN = length))
#'   # x <- within(x, vsum <- ave(get_attr(x, "validitiy"), as.list(vars), FUN = sum))
#'   x <- ave(x[[get_attr(x, "coverage")]], x[, vars, drop = TRUE], FUN = length)
#'
#'   x <- x[do.call(order, x[vars, drop=T], -x[[get_attr(x, "sample")]]), ]
#'
#'   x <- within(x, n <- ave(get_attr(x, "coverage"),
#'                           as.list(c(vars, get_attr(x, "evidence"))),
#'                           FUN = seq_along))
#'   return(x)
#' }

