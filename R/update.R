#' Update ic data to adjust denominators.
#'
#' Provide data checking and validation by comparing across official estimates,
#' administrative records, and survey estimates of immunisation coverage.
#' @param admin Object of class \code{ic.df} to containing administrative
#'   estimates of coverage.
#' @param denom Object of class \code{ic.df} containing official WUENIC
#'   denominators.
#' @param official Object of class \code{ic.df} containing official country
#'   estimates of coverage.
#' @param survey Object of class \code{ic.svy} containing annual survey
#'   estimates of coverage.
#' @param dtpCheck Logical. Should dose 3 and dose 1 of DTP vaccines be compared
#'   for each group and time? If dose 3 coverage exceeds dose 1, these records
#'   will marked \code{NA}.
#' @name ic_update
#' @export
ic_update <- function(admin, denom, official, survey, dtpCheck = TRUE){
  if(!is.ic_data(admin) || !is.ic_data(denom) ||
     !is.ic_data(official) || !is.ic_data(survey)){
    stop("Please provide valid 'ic' data.")
  }

  if(dtpCheck){
    if(all(c("DTP1", "DTP3") %in% list_vaccines(admin) &
           c("DTP1", "DTP3") %in% list_vaccines(admin) &
           c("DTP1", "DTP3") %in% list_vaccines(admin))){
      admin <- check_vd3(admin)
      official <- check_vd3(official)
      survey <- check_vd3(survey)
    }
  }

  # denominator adjustment
  a_attrs <- get_attr(admin, ic_core())
  d_attrs <- get_attr(denom, ic_core())
  o_attrs <- get_attr(official, ic_core())
  s_attrs <- get_attr(official, ic_core(survey = TRUE))

  # check for missing attributes
  if(!"dose" %in% a_attrs){ stop("Admin records missing 'dose' attribute.") }
  if(!"target" %in% d_attrs){ stop("Denom records missing 'target' attribute.") }

  # check for duplicates
  if(sum(duplicated(admin[,
                          a_attrs[c("group", "time", "vaccine")],
                          drop = TRUE])) > 0){
    stop("Duplicate records found in 'admin' dataset.")
  }

  if(sum(duplicated(official[,
                             o_attrs[c("group", "time", "vaccine")],
                             drop = TRUE])) > 0){
    stop("Duplicate records found in 'official' dataset.")
  }

  if(sum(duplicated(survey[,
                           s_attrs[c("group", "time", "vaccine")],
                           drop = TRUE])) > 0){
    stop("Duplicate records found in 'survey' dataset.")
  }

  # combine admin - denom - official estimates
  X <- merge(admin, denom,
             by.x = a_attrs[c("group", "time", "vaccine")],
             by.y = d_attrs[c("group", "time", "vaccine")],
             all.x = TRUE, attr.x = TRUE, suffixes = c("", ".d"))

  X <- merge(X, official,
             by.x = a_attrs[c("group", "time", "vaccine")],
             by.y = o_attrs[c("group", "time", "vaccine")],
             all.x = TRUE, attr.x = TRUE, suffixes = c("", ".o"))

  # recalculate with new denominator
  X[["coverage_adj"]] <- X[[a_attrs["dose"]]] / X[[paste0(d_attrs["target"], ".d")]]
  X[X$coverage_adj == 0, "coverageadj"] <- NA

  # replace values when > 100%
  # create flag
  X[["off_adj"]] <- ifelse(X[["coverage_adj"]] >= 100 &
                             !is.na(X[[paste0(o_attrs["coverage"], ".o")]]),
                           1, 0)
  X[X$off_adj == 1, "coverage_adj"] <- X[X$off_adj == 1, paste0(o_attrs["coverage"], ".o")]
  X[X$coverage_adj >= 100, "coverage_adj"] <- NA

  # clean up
  X <- X[, !grepl(".d|.y", names(X))]

  # combine and compare with Survey records
  X <- merge(X, survey,
             by.x = a_attrs[c("group", "time", "vaccine")],
             by.y = s_attrs[c("group", "time", "vaccine")],
             all.x = TRUE, attr.x = TRUE, suffixes = c("", ".s"))

  # calculate difference
  X[["diff"]] <- abs(X[[a_attrs["coverage"]]] - X[[paste0(s_attrs["coverage"], ".s")]])
  X[["coverage_adj"]] <- ifelse(X[["diff"]] > 10, X[[paste0(s_attrs["coverage"], ".s")]], X[["coverage_adj"]])

  # check again for DTP consistency vd1 > vd3
  return(X)
}


check_vd3 <- function(x){
  attrs <- get_attr(x, ic_core())
  # set DTP3 coverage to NA when > dose 1 coverage
  x <- merge(x,
             x[x[[attrs["vaccine"]]] == "DTP1", ],
             by = attrs[c("group", "time")],
             all.x = TRUE, attr.x = TRUE, suffixes = c("", ".y"))

  x[[attrs["coverage"]]] <- ifelse(x[[attrs["vaccine"]]] == "DTP3" &
                                  x[[attrs["coverage"]]] > x[[paste0(attrs["coverage"], ".y")]],
                                NA,
                                x[[attrs["coverage"]]])
  x <- x[, !grepl(".y", names(x), fixed = TRUE)]
  return(x)
}
