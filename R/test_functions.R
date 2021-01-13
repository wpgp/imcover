
as.ic_data <- function(...){
  x <- list(...)
  
  if(length(x) == 1L && inherits(x[[1L]], c("data.frame", "tbl_df"))){
    x <- x[[1L]]
  } else{
    x <- data.frame(x, stringsAsFactors = FALSE)
  }
  
  varnames <- names(x)
  if(any(is.null(varnames))){ stop("Missing data names not allowed.") }
  
  if(length(varnames) < 4){
    stop("Not enough columns. Expecting data for: group, time, vaccine, 
         (doses, target population,) and coverage.")
  }
  
  # set attributes
  class(x) <- list("ic_data", class(x))
  if(length(varnames) == 4){
    attributes(x)[c("group", "time", "vaccine", "coverage")] <- varnames[1:4]
  } else{
    attributes(x)[c("group", "time", "vaccine", 
                    "dose", "population", "coverage")] <- varnames[1:6]
  }
  
  # attributes(x)[c("group", "time", "vaccine", 
  #                 "dose", "population", "coverage")] <- varnames
  # attr(x, "group") <- varnames[1]
  # attr(x, "time") <- varnames[2]
  # attr(x, "vaccine") <- varnames[3]
  # attr(x, "doses") <- varnames[4]
  # attr(x, "population") <- varnames[5]
  # attr(x, "coverage") <- varnames[6]
  
  return(x)
}


is.ic_data <- function(object){
  return(inherits(object, "ic_data"))
}


"[.ic_data" <- function(x, i, j, ..., drop = FALSE){
  nargs <- nargs()
  attrs <- attributes(x)
  cattrs <- attrs[names(attrs) %in% ic_core()]
  # cattrs <- get_attr(x, ic_core())
  
  if(!missing(i) && nargs > 2){
    if(is.character(i)) i <- match(i, row.names(x))
  }
  
  class(x) <- setdiff(class(x), "ic_data") # drop
  if(!drop){
    coredf <- x[, unlist(cattrs)]
  }
  
  x <- if(missing(j)){
    if(nargs == 2) {
      x[i]
    }
    else x[i, , drop = drop]
  } else{
    x[i, j, drop = drop]
  }
  
  if(!drop){
    if(!missing(i) && (nargs > 2) || !missing(j)) coredf <- coredf[i,]
    
    miss_nms <- cattrs[!cattrs %in% names(x)]
    x[unlist(miss_nms)] <- coredf[, unlist(miss_nms)]
    attributes(x)[names(miss_nms)] <- unlist(miss_nms)
    
    x <- x[, c(unlist(cattrs), setdiff(names(x), unlist(cattrs)))]
    
    # remake ic data - else return unclass()
    cattrs[!cattrs %in% names(x)] <- NULL
    x <- do.call("ic_data", c(list(X=x), cattrs))
    
    # for(a in names(cattrs)){
    #   cc <- cattrs[[a]]
    #   x[[cc]] <- coredf[[cc]]
    #   attr(x, a) <- cc
    # }
    # class(x) <- list("ic_data", class(x))
  }
  
  return(x)
}


"[[<-.ic_data" <- function(x, i, value){
  attrs <- get_attr(x, ic_core(), unlist = TRUE)
  
  if(!(i %in% attrs && is.null(value))){
    x <- structure(NextMethod(), 
                   class = c("ic_data", setdiff(class(x), "ic_data")))
  } else{
    warning("Cannot drop core 'ic' data columns.", call. = FALSE, )
  }
  # x <- do.call("ic_data", c(list(X=NextMethod()), attrs))
  # x <- structure(NextMethod(), class = c("ic_data", setdiff(class(x), "ic_data")))
  
  return(x)
}


# "[[.ic_data" <- function(x, ...){
#   class(x) <- setdiff(class(x), "ic_data") # drop
#   x[[...]]
# }


"$<-.ic_data" <- function(x, i, value){
  x[[i]] <- value
}


ic_data <- function(X, group = 'iso3countrycode', time = 'year', 
                    vaccine = 'vaccine_name', coverage = 'percentcoverage',
                    dose = 'dosesadministered', population = 'targetgroup',
                    dropCols = FALSE, expand = FALSE, validate = FALSE, ...) UseMethod("ic_data")


ic_data.ic_data <- function(X, group = 'iso3countrycode', time = 'year', 
                            vaccine = 'vaccine_name', coverage = 'percentcoverage',
                            dose = 'dosesadministered', population = 'targetgroup',
                            dropCols = FALSE, expand = FALSE, validate = FALSE, ...){
  X
}

ic_data.data.frame <- function(X, 
                               group = 'iso3countrycode', time = 'year', 
                               vaccine = 'vaccine_name', coverage = 'percentcoverage',
                               dose = 'dosesadministered', population = 'targetgroup',
                               dropCols = FALSE, expand = FALSE, validate = FALSE, ...){
  
  # check data types
  if(missing(X)){
    stop("Please supply a valid dataset.")
  } 
  
  if(inherits(X, "list")){
    if(is.null(names(X))){ stop("Please supply valid data names.") }
    X <- data.frame(X, stringsAsFactors = FALSE)
  } 
  
  # check names
  if(any(lengths(list(time, vaccine, coverage, dose, population)) > 1)){
    stop("Please provide valid column names.")
  }
  
  varnames <- c(group, time, vaccine)
  chk <- !varnames %in% names(X)
  if(sum(chk) > 0){
    stop(paste(varnames[chk], collapse = " "), " not found in dataset.")
  }
  
  # get coverage
  if(!coverage %in% names(X)){
    if(dose %in% names(X) && population %in% names(X)){
      X[[coverage]] <- X[[dose]] / X[[population]] * 100
    } else{
      stop("Provide coverage data or doses and target population.") 
    }
  }
  if(!dose %in% names(X)){ dose <- NULL }
  if(!population %in% names(X)){ population <- NULL }
  
  keepnames <- c(varnames, dose, population, coverage)
  if(dropCols){
    X <- X[, keepnames]
  } else{
    othernames <- names(X)[!names(X) %in% keepnames]
    X <- X[, c(keepnames, othernames)]
  }
  
  # set attributes
  class(X) <- list("ic_data", class(X))
  attr(X, "group") <- group
  attr(X, "time") <- time
  attr(X, "vaccine") <- vaccine
  attr(X, "coverage") <- coverage
  attr(X, "dose") <- dose
  attr(X, "population") <- population
  
  if(expand){
    X <- ic_expand(X, ...)
  }
  
  if(validate){
    X <- ic_validate(X)
  }
  
  # sort
  # sort_list <- c(attr(X, "group"), attr(X, "time"), attr(X, "vaccine"))
  # X <- X[do.call(order, X[ , match(sort_list, names(X))]),]
  # row.names(X) <- seq(nrow(X))
  
  return(X)
}


ic_survey <- function(X, ..., 
                      evidence = NULL, sampleSize = NULL, 
                      minSample = 300, reduce = TRUE, biasAdjust = TRUE, 
                      dropCols = FALSE, validate = FALSE){
  if(missing(X)){
    stop("Please supply a valid dataset.")
  } else{
    X <- ic_data(X, dropCols = FALSE, expand = FALSE, validate = FALSE, ...)
  }
  
  if(!is.null(evidence)){
    
  }
  return(X)
}


ic_expand <- function(X, min = 1999, max, na.remove = TRUE){
  if(!is.ic_data(X)){
    stop("Please provide valid 'ic' data.")
  }
  # check years
  if(missing(max)){
    max <- base::max(X[[attr(X, "time")]], na.rm = TRUE)
  }
  
  if(min > max){
    stop("Invalid timespan.")
  }
  years <- min:max
  
  # expand groups
  df <- unique(X[, c(attr(X, "group"), attr(X, "vaccine"))])
  times <- rep(years, times=nrow(df))
  df <- df[rep(seq_len(nrow(df)), each=length(years)), ]
  df[[attr(X, "time")]] <- times
  # add NAs
  df <- merge(df, X, 
              by = c(attr(X, "group"), attr(X, "time"), attr(X, "vaccine")), 
              all.x = TRUE, sort = FALSE, attr.x = FALSE)
  
  # drop group x vacc where all years = NA
  if(na.remove){
    splits <- split(df, 
                    f = df[, c(attr(df, "group"), attr(df, "vaccine"))], 
                    drop = TRUE)
    empty <- lapply(splits, 
                    FUN = function(i){ all(is.na(i[[attr(df, "coverage")]])) })
    empty <- unlist(empty, use.names = FALSE)
    
    df <- do.call(rbind.data.frame, splits[which(!empty)])
  }
  
  # sort
  sort_list <- c(attr(df, "group"), attr(df, "time"), attr(df, "vaccine"))
  df <- df[do.call(order, df[ , match(sort_list, names(df))]),]
  row.names(df) <- seq(nrow(df))
  
  return(df)
}


ic_update <- function(X, compare, validate = TRUE){
  if(!is.ic_data(X)){
    stop("Please provide valid 'ic' data.")
  }
  
}


ic_filter <- function(X, yovi = 'default'){
  if(!is.ic_data(X)){
    stop("Please provide valid 'ic' data.")
  }
  
  if(missing(yovi)){
    stop("Please provide a data.frame of vaccines, or select 'default'.")
  }
  
  if(is.character(yovi)){
    if(yovi == "default"){
      yovi <- get_yovi()
    }
  } else{
    stopifnot(inherits(yovi, "data.frame"))
  }
  
  # filter
  return(X)
}


get_yovi <- function(X){
 # return full yovi table
}


ic_update_pop <- function(X, pop = "wpp", validate = TRUE){
  if(!is.ic_data(X)){
    stop("Please provide valid 'ic' data.")
  }

  TRUE
}


list_vaccines <- function(X){
  if(!inherits(X, c("ic_data", "ic"))){
    stop("Please supply valid 'ic' data")
  }
  return(unique(X[[attr(X, "vaccine")]]))
}


rbind.ic_data <- function(..., fill = TRUE, deparse.level = 1){
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
    # names(i)[names(i) %in% i_attr[[miss_nms]]] <- attrs[[miss_nms]]
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


merge.ic_data <- function(x, y, attr.x = TRUE, ...) {
  if(attr.x){
    attrs <- attributes(x)
  } else{
    attrs <- attributes(y)
  }
  
  df <- merge(as.data.frame(x), as.data.frame(y), ...)
  
  class(df) <- list("ic_data", class(df))
  
  attrs <- attrs[!names(attrs) %in% names(attributes(df))]
  if(any(!attrs %in% names(df))){ stop("Invalid columns and attributes.") }
  
  for (col in names(attrs)){ 
    attr(df, col) <- attrs[[col]] 
  }
  
  return(df)
}


add_region <- function(X, group, region = "who"){
  if(!is.ic_data(X) && missing(group)){
    stop("Please supply a grouping variable or 'ic' data.")
  }
  
  if(missing(groupVar)){
    group <- attr(X, "group")
  } else if(!is.character(group)){
    stop("Please supply grouping variable as a character name.")
  }

  if(length(group) > 1){ group <- group[1L] }
  if(!group %in% names(X)){ stop("Please supply a valid column name.")}
  
  X[["region"]] <- lookup_region(X[[group]], region)
}


lookup_region <- function(code, region){
  return(region_tbl[region_tbl$code == code,][[region]])
}


ic_core <- function(){ return(c("group","time","vaccine","coverage","dose","population")) }


get_attr <- function(X, attrs, unlist = TRUE){
  if(missing(attrs)) stop("Must provide one or more attribute names to match")
  
  if(unlist){
    return(unlist(attributes(X)[which(names(attributes(X)) %in% attrs)]))
  } else{
    return(attributes(X)[which(names(attributes(X)) %in% attrs)])
  }
}
