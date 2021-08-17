#' Expand ic data to include all potential years.
#'
#' Expand the time dimension of an ic dataset to include the full set of
#' potential observations.
#' @param X Object of class \code{ic.df} to expand.
#' @param min Start of the potential time period. Default is 1999.
#' @param max End of the potential time period.
#' @param na.remove Should country-vaccine combinations which have no coverage
#'   data be removed? Default is \code{TRUE}
#' @details If \code{max} is missing, the maximum time value observed in the
#'   dataset will be used.
#' @return An object of type \code{ic.df} which has additional rows of \code{NA}
#'   values added to expand the dataset to min and max time points.
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

  # drop other years
  X <- X[X[[attr(X, "time")]] >= min & X[[attr(X, "time")]] <= max, ]

  # expand groups
  df <- unique(X[, c(attr(X, "country"), attr(X, "vaccine")), drop = TRUE])
  times <- rep(years, times=nrow(df))
  df <- df[rep(seq_len(nrow(df)), each=length(years)), ]
  df[[attr(X, "time")]] <- times
  # add NAs
  df <- merge(X, df, # merge.ic.df
              by = c(attr(X, "country"), attr(X, "time"), attr(X, "vaccine")),
              all.y = TRUE, sort = FALSE, attr.x = TRUE)

  # drop country x vacc where all years = NA
  if(na.remove){
    splits <- split(df,
                    f = df[, c(attr(X, "country"), attr(X, "vaccine")), drop = TRUE],
                    drop = TRUE)
    empty <- lapply(splits,
                    FUN = function(i){ all(is.na(i[[attr(X, "coverage")]])) })
    empty <- unlist(empty, use.names = FALSE)

    df <- do.call(rbind.data.frame, splits[which(!empty)])
  }

  # sort
  sort_list <- c(attr(X, "country"), attr(X, "time"), attr(X, "vaccine"))
  df <- df[do.call(order, df[ , match(sort_list, names(df))]),]
  row.names(df) <- seq(nrow(df))

  # df <- do.call("ic_data", c(list(X=df), attributes(X)))
  return(df)
}
