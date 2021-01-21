

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
    miss_nms <- names(which(attrs != i_attr))

    for(n in miss_nms){
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


merge.ic.df <- function(x, y, attr.x = TRUE, ...) {
  if(attr.x){
    attrs <- attributes(x)
  } else{
    attrs <- attributes(y)
  }

  df <- merge(as.data.frame(x), as.data.frame(y), ...)

  class(df) <- list("ic.df", class(df))

  attrs <- attrs[!names(attrs) %in% names(attributes(df))]
  if(any(!attrs %in% names(df))){ stop("Invalid columns and attributes.") }

  for (col in names(attrs)){
    attr(df, col) <- attrs[[col]]
  }

  return(df)
}
