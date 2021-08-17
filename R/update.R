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
#' @param dtpCheck Should dose 3 and dose 1 of DTP vaccines be compared for each
#'   country and time? If dose 3 coverage exceeds dose 1, these records will
#'   marked \code{NA}. Default is \code{TRUE}.
#' @param limit Maximum acceptable difference between survey and administrative
#'   coverage records. When the difference exceeds the limit, the survey
#'   estimate is preferred. Default is 10 percent.
#' @seealso \code{\link[imcover]{ic_data}}
#' @name ic_update
#' @export
ic_update <- function(admin, denom, official, survey, dtpCheck = TRUE, limit = 10){
  if(!is.ic_data(admin) || !is.ic_data(denom) ||
     !is.ic_data(official) || !is.ic_data(survey)){
    stop("Please provide valid 'ic' data.")
  }

  if(dtpCheck){
    if(!all(c("DTP1", "DTP3") %in% list_vaccines(admin) &
           c("DTP1", "DTP3") %in% list_vaccines(official) &
           c("DTP1", "DTP3") %in% list_vaccines(survey))){
      stop("DTP vaccines not present in data.") }

      admin <- check_vd3(admin)
      official <- check_vd3(official)
      survey <- check_vd3(survey)
  }

  # denominator adjustment
  a_attrs <- get_attr(admin, ic_core())
  d_attrs <- get_attr(denom, ic_core())
  o_attrs <- get_attr(official, ic_core())
  s_attrs <- get_attr(survey, ic_core(survey = TRUE))

  # check for missing attributes
  if(!"dose" %in% names(a_attrs)){ stop("Admin records missing 'dose' attribute.") }
  if(!"population" %in% names(d_attrs)){ stop("Denom records missing 'population' attribute.") }

  # check for duplicates
  if(sum(duplicated(admin[,
                          a_attrs[c("country", "time", "vaccine")],
                          drop = TRUE])) > 0){
    stop("Duplicate records found in 'admin' dataset.")
  }

  if(sum(duplicated(official[,
                             o_attrs[c("country", "time", "vaccine")],
                             drop = TRUE])) > 0){
    stop("Duplicate records found in 'official' dataset.")
  }

  if(sum(duplicated(survey[,
                           s_attrs[c("country", "time", "vaccine")],
                           drop = TRUE])) > 0){
    stop("Duplicate records found in 'survey' dataset.")
  }

  # combine admin - denom - official estimates
  names(denom) <- paste0(names(denom), ".d")
  X <- merge(admin, denom,
             by.x = a_attrs[c("country", "time", "vaccine")],
             by.y = paste0(d_attrs[c("country", "time", "vaccine")], ".d"),
             all.x = TRUE, attr.x = TRUE)

  names(official) <- paste0(names(official), ".o")
  X <- merge(X, official,
             by.x = a_attrs[c("country", "time", "vaccine")],
             by.y = paste0(o_attrs[c("country", "time", "vaccine")], ".o"),
             all.x = TRUE, attr.x = TRUE)

  # recalculate with new denominator
  X[["coverage_adj"]] <- X[[a_attrs["dose"]]] / X[[paste0(d_attrs["population"], ".d")]] * 100
  X$coverage_adj <- ifelse(X$coverage_adj == 0, NA, X$coverage_adj)

  # replace values when new coverage > 100%
  # create flag
  X[["off_adj"]] <- ifelse(X[["coverage_adj"]] >= 100 &
                             !is.na(X[[paste0(o_attrs["coverage"], ".o")]]),
                           1, 0)

  # combine multi-dose checks/adjustments
  # how to handle NA values?
  if(dtpCheck){
    mdadj <- X[X[[a_attrs['vaccine']]] %in% c("DTP1", "DTP3"),
               c(a_attrs[c("country", "time", "vaccine")], "off_adj"), drop = TRUE]
    i <- which(X[[a_attrs['vaccine']]] %in% c("DTP1", "DTP3"))
    mdadj <- ave(mdadj$off_adj, mdadj[a_attrs[c("country", "time")]], FUN = max)
    X[i, "off_adj"] <- mdadj
  }
  # adjustment
  X$coverage_adj <- ifelse(X$off_adj == 1, X[[paste0(o_attrs["coverage"], ".o")]], X$coverage_adj)
  X[X$coverage_adj >= 100 & !is.na(X$coverage_adj), "coverage_adj"] <- NA

  # clean up
  X <- X[, !grepl("\\.d|\\.o", names(X))]

  ###
  # combine and compare with Survey records
  names(survey) <- paste0(names(survey), ".s")
  X <- merge(X, survey,
             by.x = a_attrs[c("country", "time", "vaccine")],
             by.y = paste0(s_attrs[c("country", "time", "vaccine")], ".s"),
             all.x = TRUE, attr.x = TRUE)

  # calculate difference
  X[["diff"]] <- abs(X[[a_attrs["coverage"]]] - X[[paste0(s_attrs["coverage"], ".s")]])
  # create flag
  X[["surv_adj"]] <- ifelse(abs(X[["diff"]]) > limit &
                              !is.na(X[[paste0(s_attrs["coverage"], ".s")]]),
                            1, 0)
  # handle multi-dose comparisons
  if(dtpCheck){
    mdadj <- X[X[[a_attrs['vaccine']]] %in% c("DTP1", "DTP3"),
               c(a_attrs[c("country", "time", "vaccine")], "surv_adj"), drop = TRUE]
    i <- which(X[[a_attrs['vaccine']]] %in% c("DTP1", "DTP3"))
    mdadj <- ave(mdadj$surv_adj, mdadj[a_attrs[c("country", "time")]], FUN = max)
    X[i, "surv_adj"] <- mdadj
  }
  # adjustment
  X$coverage_adj <- ifelse(X$surv_adj == 1, X[[paste0(s_attrs["coverage"], ".s")]], X$coverage_adj)

  # clean up
  X <- X[, !grepl("\\.s", names(X))]

  # check again for DTP consistency vd1 > vd3
  X <- check_vd3(X)

  # set coverage attribute
  attributes(X)["coverage"] <- "coverage_adj"
  return(X)
}


check_vd3 <- function(x){
  attrs <- get_attr(x, ic_core())
  # set DTP3 coverage to NA when > dose 1 coverage
  x <- merge(x,
             x[x[[attrs["vaccine"]]] == "DTP1", ],
             by = attrs[c("country", "time")],
             all.x = TRUE, attr.x = TRUE, suffixes = c("", ".y"))

  x[[attrs["coverage"]]] <- ifelse(x[[attrs["vaccine"]]] == "DTP3" &
                                  x[[attrs["coverage"]]] > x[[paste0(attrs["coverage"], ".y")]],
                                NA,
                                x[[attrs["coverage"]]])
  x <- x[, !grepl(".y", names(x), fixed = TRUE)]
  return(x)
}
