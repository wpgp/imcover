#' Combine ic data objects by rows
#'
#' Take a sequence of \code{ic.df} data-frame arguments and combine by rows.
#' @param ... Objects to bind together.
#' @param fill Should missing columns be added to data objects so that they can
#'   be combined? Default is \code{TRUE}.
#' @param deparse.level See \code{\link{rbind}} Default is 1.
#' @details The attributes of the ic data are used to match columns. Column
#'   names are taken from the first object.
#' @return An object of type \code{ic.df}.
#' @name rbind
#' @export
rbind.ic.df <- function(..., fill = TRUE, deparse.level = 1){
  dots <- list(...)
  dots <- dots[!sapply(dots, is.null)]

  if(length(dots) > 1L){
    is_ic <- vapply(dots[-1L], function(i) is.ic_data(i), TRUE)
    if(!all(is_ic)) stop("All objects must be 'ic' data.")
  }

  # match/rename core attributes + columns
  attrs <- get_attr(dots[[1L]], ic_core())
  dots <- lapply(dots, function(i){
    i_attr <- get_attr(i, ic_core())
    # miss_nms <- names(which(attrs != i_attr))
    #
    # for(n in miss_nms){
    #   names(i)[names(i) == i_attr[[n]]] <- attrs[[n]]
    #   attr(i, n) <- attrs[[n]]
    # }

    match_attr <- intersect(names(attrs), names(i_attr))
    tofix <- which(!get_attr(i, match_attr) %in% attrs[match_attr])

    for(n in match_attr[tofix]){
      names(i)[names(i) == i_attr[[n]]] <- attrs[[n]]
      attr(i, n) <- attrs[[n]]
    }

    return(i)
  })

  if(fill){
    nms <- unique(unlist(lapply(dots, function(i) names(i)), use.names = FALSE))
    dots <- lapply(dots,
                   function(i){ i[, setdiff(nms, names(i))] <- NA; return(i) })
  }

  df <- do.call("rbind.data.frame", dots)
  return(df)
}


#' Merge two ic data objects
#' Merge two \code{ic.df} data objects by common columns while preserving the
#' attributes
#' @param x,y \code{ic.df} objects to merge
#' @param attr.x Logical. Should the attributes (i.e. core field names) be based
#'   on the ic object 'x' or 'y' (when \code{FALSE})?
#' @param ... arguments to be based on \code{merge}.
#' @export
merge.ic.df <- function(x, y, attr.x = TRUE, ...) {
  if(attr.x){
    attrs <- attributes(x)
    drops <- attributes(y)
  } else{
    attrs <- attributes(y)
    drops <- attributes(x)
  }
  attrs <- attrs[names(attrs) %in% ic_core(survey = TRUE)]
  drops <- drops[names(drops) %in% ic_core(survey = TRUE)]
  drops <- setdiff(unlist(drops), unlist(attrs))

  df <- merge(as.data.frame(x), as.data.frame(y), ...)
  df <- df[,!names(df) %in% drops]

  if(any(!unlist(attrs) %in% names(df))){
    stop("Invalid columns and attributes.")
  }
  # attributes(df)[names(attrs)] <- attrs
  # class(df) <- list("ic.df", class(df))

  df <- do.call("ic_data", c(list(X=df), attrs))
  stopifnot(is.ic_data(df))
  # print(attributes(df))

  return(df)
}

