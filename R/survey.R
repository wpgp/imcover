

#' Create an ic data object
#' @description Create an ic data object from survey data on immunisation
#'   coverage.
#' @param X Any object with immunisation coverage survey data to convert into an
#'   ic object.
#' @param ... Additional parameters for names of core data elements to pass on
#'   \code{ic_data}. See details and \code{ic_data}.
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

  # correct dose 3 recall
  if(biasAdjust){
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

    if(!is.logical(X[[validity]])){
      stop("Validity indicator must be a boolean (logical) value.")
    }
    # potentially missing/malformed
    X[[evidence]] <- trimws(tolower(as.character(X[[evidence]])))

    # identify groupings
    X <- mark_survey(X)

    # process data by vaccine
    dd <- split(X, f = data[[get_attr(X, "vaccine")]])

    for(v in c("DTP", "PCV")){ # change/expand here
      vd3 <- dd[[paste0(v, 3)]]
      vd1 <- dd[[paste0(v, 1)]]
      # mark for keepers
      vd3 <- mark_survey(vd3)
      vd1 <- mark_survey(vd1)

      # split data to process
      # vd3_coh <- vd3[vd3[[get_attr(vd3, "evidence")]] == "card or history" &
      #                  !is.na(vd3[[get_attr(vd3, "evidence")]]), ]
      #
      # vd3_c <- vd3[vd3[[get_attr(vd3, "evidence")]] == "card" &
      #                !is.na(vd3[[get_attr(vd3, "evidence")]]), ]

      vd1_coh <- vd1[vd1[[get_attr(vd1, "evidence")]] == "card or history" &
                       !is.na(vd1[[get_attr(vd1, "evidence")]]), ]

      vd1_c <- vd1[vd1[[get_attr(vd1, "evidence")]] == "card" &
                     !is.na(vd1[[get_attr(vd1, "evidence")]]), ]
      # merge testing
      # c1 <- dtp1[dtp1$evidence == "Card" & (dtp1$vv | dtp1$N == 1) & dtp1$n == 1,]
      # ch1 <- dtp1[dtp1$evidence == "Card or History" & (dtp1$vv | dtp1$N ==1) & dtp1$n == 1,]
      #
      # dtp1 <- merge(c1, ch1, by = c("ISO3", "surveyNameEnglish", "cohortYear"))
#
#       vd_update <- merge(vd3_coh[(vd3_coh[[get_attr(vd3_coh, "valid")]] | vd3_coh$N == 1) & vd3_coh$n == 1, ],
#                          vd3_c[(vd3_c[[get_attr(vd3_c, "valid")]] | vd3_c$N == 1) & vd3_c$n == 1, ],
#                          by = c(get_attr(vd3_coh, c("group", "survey", "time"))))

      vd1_update <- merge(vd1_coh[(vd1_coh[[get_attr(vd1_coh, "valid")]] | vd1_coh$N == 1) & vd1_coh$n == 1, ],
                          vd1_c[(vd1_c[[get_attr(vd1_c, "valid")]] | vd1_c$N == 1) & vd1_c$n == 1, ],
                          by = c(get_attr(vd1_coh, c("group", "survey", "time"))))
      vd1_update$adj_factor <- vd1_update[[paste0(get_attr(vd1_coh, "coverage"), ".x")]] * vd1_update[[paste0(get_attr(vd1_c, "coverage"), ".y")]]
    }

    # DTP
    dtp3 <- dd[["DTP3"]]
    dtp1 <- dd[["DTP1"]]

  }

  # drop group markings

  return(X)
}


mark_survey <- function(x){
  vars <- get_attr(x, c("group", "survey", "time", "vaccine"))
  x <- within(x, N <- ave(get_attr(x, "coverage"), as.list(vars), FUN = length))
  # x <- within(x, vsum <- ave(get_attr(x, "validitiy"), as.list(vars), FUN = sum))

  x <- x[do.call(order, x[vars, drop=T], -x[[get_attr(x, "sample")]]), ]

  x <- within(x, n <- ave(get_attr(x, "coverage"),
                          as.list(c(vars, get_attr(x, "evidence"))),
                          FUN = seq_along))
  return(x)
}
