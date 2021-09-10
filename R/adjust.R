
#' Adjust 'ic' coverage data
#'
#' Make adjustments to vaccination coverage to ensure reported coverages
#' preserve the ratio of multi-dose vaccines and that coverage is >0% and <100%.
#' @param X Object of class \code{ic.df} to adjust
#' @param ratio_adj Logical. Should multi-dose vaccine coverages be adjusted by
#'   their ratio when coverage > 100%? Default is \code{FALSE}.
#' @param numerator Character vector of vaccine identifiers (with dose numbers)
#'   to calculate the ratio adjustment. Default is 'DTP3'.
#' @param denominator Character vector of vaccine identifiers (with dose
#'   numbers) to calculate the ratio adjustment. Default is 'DTP1'.
#' @param coverage_adj Logical. Should coverage percentage be changed to be >0%
#'   and <100%? Default is \code{TRUE}.
#' @details This function provides minimal adjustments to immunisation data to
#'   ensure that the coverages are suitable for statistical modelling. Ratio
#'   adjustments help to ensure that multi-dose vaccines report a coverage that
#'   is less-than-or-equal-to the coverage of the first dose.
#'
#'   Coverage adjustment accounts for inaccurate target populations that result
#'   is coverages of >100%. By maintaining a coverage between 0 and 100% a
#'   logistic transformation can be used in the statistical model.
#' @seealso \code{\link[imcover]{ic_data}}
#' @name ic_adjust
#' @export
ic_adjust <- function(X, ratio_adj = FALSE,
                      numerator = 'DTP3', denominator = 'DTP1',
                      coverage_adj = TRUE){
  if(missing(X) || !inherits(X, 'ic.df')){
    stop('Please provide valid `ic` data')
  }

  if(any(!(c(numerator, denominator) %in% list_vaccines(X)))){
    stop('Vaccines to adjust not found in ic data.')
  }

  if(length(numerator) != length(denominator)){
    stop('Length of numerator and denominator for ratio adjustment do not match.')
  }

  if(ratio_adj){
    for(i in 1:length(numerator)){
      num <- as.data.frame(X[X[[attr(X, 'vaccine')]] == numerator[[i]], ])
      den <- as.data.frame(X[X[[attr(X, 'vaccine')]] == denominator[[i]], ])
      # drop vacc records
      X <- X[X[[attr(X, 'vaccine')]] != numerator[[i]], ]

      df <- merge(num,
                  den[, get_attr(X, c('country', 'time', 'source', 'coverage'))],
                  by = get_attr(X, c('country', 'time', 'source')))

      df$above100 <- ifelse(df[[paste0(attr(X, 'coverage'), '.y')]] > 100,
                            1, 0)
      # calculate the ratio
      df$d_ratio <- ifelse(df$above100 == 1,
                           df[[paste0(get_attr(X, 'coverage'), '.x')]] / df[[paste0(get_attr(X, 'coverage'), '.y')]],
                           NA) # dose 3 / dose 1

      # adjustment of dose 1 values
      df[[paste0(get_attr(X, 'coverage'), '.y')]] <- ifelse(df[[paste0(get_attr(X, 'coverage'), '.y')]] >= 100,
                                                              99.9,
                                                              df[[paste0(get_attr(X, 'coverage'), '.y')]])

      df[[paste0(get_attr(X, 'coverage'), '.y')]] <- ifelse(df[[paste0(get_attr(X, 'coverage'), '.y')]] < 1,
                                                              0.1,
                                                              df[[paste0(get_attr(X, 'coverage'), '.y')]])

      # correction of dose 3 values using ratio and rounded dose 1 coverage
      df$coverage.x <- ifelse(!is.na(df$d_ratio),
                               df$d_ratio * df$coverage.y,
                               df$coverage.y)

      # reassemble dataset
      vd3 <- df[, c(intersect(names(df), names(X)),
                    names(df)[grepl(".x", names(df), fixed = T)])]
      names(vd3) <- gsub(".x", "", names(vd3), fixed = T)
      # convert to ic
      vd3 <- do.call("ic_data", c(list(X=vd3), get_attr(X, ic_core())))

      X <- rbind(X, vd3)
    }
  }

  if(coverage_adj){
    X[[attr(X, 'coverage')]] <- ifelse(X[[attr(X, 'coverage')]] >= 100,
                                       99.9,
                                       X[[attr(X, 'coverage')]])

    X[[attr(X, 'coverage')]] <- ifelse(X[[attr(X, 'coverage')]] < 1,
                                       0.1,
                                       X[[attr(X, 'coverage')]])
  }

  return(X)
}


#' Calculate the ratio between multi-dose vaccines
#'
#' @param X Object of class \code{ic.df} to adjust
#' @param numerator Character vector of vaccine identifiers (with dose numbers)
#'   to calculate the ratio adjustment. Default is 'DTP3'.
#' @param denominator Character vector of vaccine identifiers (with dose
#'   numbers) to calculate the ratio adjustment. Default is 'DTP1'.
#' @name ic_adjust
#' @export
ic_ratio <- function(X, numerator = 'DTP3', denominator = 'DTP1'){

  if(missing(X) || !inherits(X, 'ic.df')){
    stop('Please provide valid `ic` data')
  }

  if(any(!(c(numerator, denominator) %in% list_vaccines(X)))){
    stop('Vaccines to adjust not found in ic data.')
  }

  if(length(numerator) != length(denominator)){
    stop('Length of numerator and denominator for ratio adjustment do not match.')
  }

  for(i in 1:length(numerator)){
    num <- as.data.frame(X[X[[attr(X, 'vaccine')]] == numerator[[i]], ])
    den <- as.data.frame(X[X[[attr(X, 'vaccine')]] == denominator[[i]], ])
    # drop vacc records
    X <- X[X[[attr(X, 'vaccine')]] != numerator[[i]], ]

    df <- merge(num,
                den[, get_attr(X, c('country', 'time', 'source', 'coverage'))],
                by = get_attr(X, c('country', 'time', 'source')))

    # calculate the ratio - multiply by 100 to match coverage percentages
    df$d_ratio <- df[[paste0(get_attr(X, 'coverage'), '.x')]] / df[[paste0(get_attr(X, 'coverage'), '.y')]] * 100

    df$d_ratio <- ifelse(is.na(df$d_ratio) | df$d_ratio < 1, 0.1, df$d_ratio)
    df$d_ratio <- ifelse(df$d_ratio >= 100, 99.9, df$d_ratio)
    df$coverage <- df$d_ratio

    # reassemble dataset
    vd3 <- df[, c(intersect(names(df), names(X)))]
    # convert to ic
    vd3 <- do.call("ic_data", c(list(X=vd3), get_attr(X, ic_core())))

    X <- rbind(X, vd3)
  }

  # add non-core attribute to flag data have been modified
  attr(X, 'numerator') <- numerator
  attr(X, 'denominator') <- denominator

  return(X)
}
