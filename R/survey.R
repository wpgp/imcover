

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
ic_survey <- function(X, ..., survey = "surveyname",
                      evidence = NULL, sampleSize = NULL,
                      minSample = NULL, reduce = TRUE, biasAdjust = TRUE,
                      dropCols = FALSE, validate = FALSE){
  UseMethod('ic_survey')
}


#' @name ic_survey
#' @export
ic_survey.ic.df <- function(X, ...,
                            survey = "survey_name", sample = "sample_size",
                            evidence = "evidence", validity = "validity",
                            minSample = NULL, biasAdjust = TRUE,
                            dropCols = FALSE){

  X # add processing steps for bias and sample corrections
}


#' @name ic_survey
#' @export
ic_survey.data.frame <- function(X, ...,
                                 survey = "survey_name", sample = "sample_size",
                                 evidence = "evidence", validity = "validity",
                                 minSample = NULL, biasAdjust = TRUE,
                                 dropCols = FALSE){
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
  # potentially missing/malformed
  X[[evidence]] <- trimws(tolower(as.character(X[[evidence]])))
  X[[validity]] <- trimws(tolower(as.character(X[[validity]])))
  # sort
  attrs <- get_attr(X, attrs = ic_core(survey = TRUE), unlist = FALSE)
  corenames <- unlist(attrs)
  X <- X[, c(corenames, setdiff(names(X), corenames))]

  # correct dose 3 recall
  if(biasAdjust){
    # check for missing validity info? some empty strings. drop?
    if(any(!X[[validity]] %in% c("crude", "valid"))){
      stop(paste0("Validity must be 'crude' or 'valid'. Found: ",
                  paste(unique(X[[validity]]), collapse = ", ")),
           call. = FALSE)
    }

    X[["adj_factor"]] <- X[["coverage_adj"]] <- NA

    for(v in c("DTP", "PCV")){
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
      # vd3.coh <- merge(vd3[vd3[[attrs$evidence]] == "card or history", corenames],
      #                  vd_adj[, c(corenames[c('group','time','survey')],
      #                             c("adj_factor","coverage_adj"))],
      #                  by = corenames[c('group','time','survey')],
      #                  all.x = TRUE)
      # vd3.coh[[attrs$coverage]] <- vd3.coh[["coverage_adj"]]
      # add rows back
      # vd1[["adj_factor"]] <- vd3[["adj_factor"]] <- NA
      # vd1[["coverage_adj"]] <- vd3[["coverage_adj"]] <- NA
      # vd <- rbind(vd1, vd3.coh, vd3[vd3[[attrs$evidence]] != "card or history",
      #                               c(corenames, "adj_factor", "coverage_adj")])
      # dropping v and dose == 3
      # X <- rbind(X, vd3)
    }






#
#     # subset vaccines to process
#     dtp1 <- as.data.frame(X[X[[corenames['vaccine']]] == "DTP1", ])
#     dtp3 <- as.data.frame(X[X[[get_attr(X, "vaccine")]] == "DTP3", ])
#     # get attributes
#     attrs <- attributes(X)
#
#     # reshape to wide
#     dtp1 <- reshape(dtp1, direction = "wide",
#                     idvar = get_attr(X, c("group","survey","time",
#                                           "vaccine","evidence")),
#                     timevar = get_attr(X, "validity"))
#     # find preferred data
#     timevar <- paste(attrs$coverage, "valid", sep = ".")
#     vd1[[attrs$validity]] <- ifelse(is.na(vd1[[timevar]]),
#                                      "crude", "valid")
#
#     vd1[[attrs$coverage]] <- ifelse(is.na(vd1[[timevar]]),
#                                      vd1[[paste(attrs$coverage, "crude")]],
#                                      vd1[[paste(attrs$coverage, "valid")]])
#
#     vd1[[attrs$sample]] <- ifelse(is.na(vd1[[timevar]]),
#                                      vd1[[paste(attrs$sample, "crude", sep = ".")]],
#                                      vd1[[paste(attrs$sample, "valid", sep = ".")]])
#
#     # merge c and coh
#     vd1 <- merge(dtp1[dtp1[[attrs$evidence]] == "card or history", corenames],
#                  dtp1[dtp1[[attrs$evidence]] == "card", corenames],
#                  by = corenames[c('group','time','survey')],
#                  all.x = TRUE, suffixes = c(".coh", ".c"))
#
#
#
#
#
#     dtp3 <- reshape(dtp3, direction = "wide",
#                     idvar = get_attr(X, c("group","survey","time",
#                                           "vaccine","evidence")),
#                     timevar = get_attr(X, "validity"))
    # reshape again? to account for evidence + validity
    # process columns to find preferred data

    # convert to ic data
    # either remerge with other vacc data or update full dataset (how to ID?)


#     ### bias adjustment
#     redf <- reshape(as.data.frame(data),
#             direction = "wide",
#             idvar = c("iso3","surveynameenglish","vaccine", "evidence"),
#             timevar = "validity")
#     # find preferred
#
#     reshape(redf,
#             direction = "long",
#             varying = paste(c("cohortyear","coverage","sample_size","x","collectbegin","collectend","cardsseen","agevaccination","ageinterview","sex"),
#                             rep(c("crude","valid"), times=10), sep="."))
#
#
#     # identify groupings
#     X <- mark_survey(X)
#
#     # process data by vaccine
#     dd <- split(X, f = data[[get_attr(X, "vaccine")]])
#
#     for(v in c("DTP", "PCV")){ # change/expand here
#       vd3 <- dd[[paste0(v, 3)]]
#       vd1 <- dd[[paste0(v, 1)]]
#       # mark for keepers
#       vd3 <- mark_survey(vd3)
#       vd1 <- mark_survey(vd1)
#
#       # split data to process
#       # vd3_coh <- vd3[vd3[[get_attr(vd3, "evidence")]] == "card or history" &
#       #                  !is.na(vd3[[get_attr(vd3, "evidence")]]), ]
#       #
#       # vd3_c <- vd3[vd3[[get_attr(vd3, "evidence")]] == "card" &
#       #                !is.na(vd3[[get_attr(vd3, "evidence")]]), ]
#
#       vd1_coh <- vd1[vd1[[get_attr(vd1, "evidence")]] == "card or history" &
#                        !is.na(vd1[[get_attr(vd1, "evidence")]]), ]
#
#       vd1_c <- vd1[vd1[[get_attr(vd1, "evidence")]] == "card" &
#                      !is.na(vd1[[get_attr(vd1, "evidence")]]), ]
#       # merge testing
#       # c1 <- dtp1[dtp1$evidence == "Card" & (dtp1$vv | dtp1$N == 1) & dtp1$n == 1,]
#       # ch1 <- dtp1[dtp1$evidence == "Card or History" & (dtp1$vv | dtp1$N ==1) & dtp1$n == 1,]
#       #
#       # dtp1 <- merge(c1, ch1, by = c("ISO3", "surveyNameEnglish", "cohortYear"))
# #
# #       vd_update <- merge(vd3_coh[(vd3_coh[[get_attr(vd3_coh, "valid")]] | vd3_coh$N == 1) & vd3_coh$n == 1, ],
# #                          vd3_c[(vd3_c[[get_attr(vd3_c, "valid")]] | vd3_c$N == 1) & vd3_c$n == 1, ],
# #                          by = c(get_attr(vd3_coh, c("group", "survey", "time"))))
#
#       vd1_update <- merge(vd1_coh[(vd1_coh[[get_attr(vd1_coh, "valid")]] | vd1_coh$N == 1) & vd1_coh$n == 1, ],
#                           vd1_c[(vd1_c[[get_attr(vd1_c, "valid")]] | vd1_c$N == 1) & vd1_c$n == 1, ],
#                           by = c(get_attr(vd1_coh, c("group", "survey", "time"))))
#       vd1_update$adj_factor <- vd1_update[[paste0(get_attr(vd1_coh, "coverage"), ".x")]] * vd1_update[[paste0(get_attr(vd1_c, "coverage"), ".y")]]
#     }
#
#     # DTP
#     dtp3 <- dd[["DTP3"]]
#     dtp1 <- dd[["DTP1"]]

  }

  # drop group markings

  return(X)
}

#' Find survey groups
mark_survey <- function(x){
  vars <- get_attr(x, c("group", "survey", "time", "vaccine"))
  # x <- within(x, N <- ave(get_attr(x, "coverage"), as.list(vars), FUN = length))
  # x <- within(x, vsum <- ave(get_attr(x, "validitiy"), as.list(vars), FUN = sum))
  x <- ave(x[[get_attr(x, "coverage")]], x[, vars, drop = TRUE], FUN = length)

  x <- x[do.call(order, x[vars, drop=T], -x[[get_attr(x, "sample")]]), ]

  x <- within(x, n <- ave(get_attr(x, "coverage"),
                          as.list(c(vars, get_attr(x, "evidence"))),
                          FUN = seq_along))
  return(x)
}

