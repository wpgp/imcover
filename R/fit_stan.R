
#' Multi-source immunisation coverage model with Stan
#'
#' @param X Object of \code{ic.df} for analysis
#' @param region Logical. Should region-specific models be generated? Default
#'   is \code{TRUE}.
#' @param verbose Logical. Should messages be displayed? Default is \code{TRUE}.
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'
#' @name ic_fit
#' @export
ic_fit <- function(X, region = TRUE, verbose = TRUE, ...){
  UseMethod("ic_fit")
}


#' @name ic_fit
#' @export
ic_fit.ic.df <- function(X, region = TRUE, verbose = TRUE, ...){

  # check data
  X <- X[!is.na(X[[get_attr(X, 'coverage')]]), ]

  if(region){
    X <- X[!is.na(X[[get_attr(X, 'region')]]), ]

    # get regions
    regions <- sort(unique(X[[get_attr(X, 'region')]]))
    regions <- regions[!is.na(regions)]
  }

  # calculation
  if(!region || length(regions) == 1){
    # calculate
    out <- multi_lik_stan(X, verbose, ...)
  } else{
    # main processing loop
    out <- lapply(regions, function(r){
      if(verbose) print(r)
      # subset
      dat <- X[X[[get_attr(X, 'region')]] == r, ]
      # calculate
      fit <- multi_lik_stan(dat, verbose, ...)

      return(fit)
    })
    names(out) <- regions
    class(out) <- list("iclist", class(out))
  }

  return(out)
}


#' Multi-source immunisation coverage model with Stan
#'
#' @param X Object of \code{ic.df} for analysis
#' @param verbose Logical. Should messages be displayed? Default is \code{TRUE}.
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'
#' @name multi_lik_stan
#' @keywords internal
multi_lik_stan <- function(X, verbose = TRUE, ...) {
  if(!is.ic_data(X)) stop("Please provide valid 'ic' data.")

  # data prep
  vax_data <- ic_to_stan(X, make_index = TRUE)

  # call stan model
  out <- rstan::sampling(stanmodels$multi_lik,
                         data = vax_data,
                         show_messages = verbose,
                         ...)

  # calculate prediction ('mu')
  posterior <- t(as.data.frame(out, 'mu'))
  posterior <- invlogit(posterior)

  # create an index to the 'mu' parameter
  mu_idx <- strsplit(gsub('\\[|\\]', '',
                          regmatches(row.names(posterior),
                                     gregexpr("\\[.*?\\]",
                                              row.names(posterior)))),
                     ',', fixed = T)
  mu_idx <- do.call(rbind.data.frame, mu_idx)
  names(mu_idx) <- c("country", "time", "vaccine")
  mu_idx <- as.data.frame(lapply(mu_idx, function(x) as.numeric(x)))

  # name the indices
  mu_names <- cbind.data.frame(
    'country' = levels(factor(X[[get_attr(X, 'country')]]))[mu_idx$country],
    'time' = as.numeric(levels(factor(X[[get_attr(X, 'time')]]))[mu_idx$time]),
    'vaccine' = sort(unique(X[[get_attr(X, 'vaccine')]]))[mu_idx$vaccine]
  )

  # ratio adjustment
  if(!is.null(attr(X, 'numerator'))){
    numerator <- attr(X, 'numerator')
    denominator <- attr(X, 'denominator')

    # apply ratio to dose 1 to calc new dose 3
    for(i in 1:length(numerator)){
      num <- posterior[which(mu_names$vaccine == numerator[i]), ]
      den <- posterior[which(mu_names$vaccine == denominator[i]), ]

      num <- den * num
      # reassemble
      posterior[which(mu_names$vaccine == numerator[i]), ] <- num
    }
  }

  # convert to coverage percentage
  posterior <- posterior * 100

  # apply labels
  posterior <- cbind(mu_names, posterior)

  out <- list('fit' = out, 'posterior' = posterior, 'data' = vax_data)
  class(out) <- list('icfit', class(out))

  return(out)
}


#' Prepare 'ic' data for Stan model
#'
#' @param X Object of \code{ic.df} to be converted
#' @param make_index Logical. Should indices for e.g. 'country' be created?
#' @return A list with named objects suitable for use as Stan data.
#'   Default is \code{TRUE}.
#' @keywords internal
ic_to_stan <- function(X, make_index = TRUE){
  stopifnot(is.ic_data(X))

  # swap names
  X <- swap_names(X)
  X <- data.frame(X)

  # sort
  X <- X[order(X$source, X$country, X$vaccine, X$time), ]

  # modify outcome
  X$coverage <- X$coverage / 100
  X$cov.logit <-  logit(X$coverage)  # log(X$coverage / (1 - X$coverage))

  # create index values
  X$v <- as.numeric(factor(X$vaccine))
  X$i <- as.numeric(factor(X$country))
  X$t <- as.numeric(factor(X$time)) # temporal trend per country
  X$s <- as.numeric(factor(X$source))

  vax_dat <- list(y = X$cov.logit,
                  i = X$i,
                  j = X$v,
                  t = X$t,
                  N = nrow(X),
                  N_a = sum(tolower(X$source) == "admin"),
                  N_o = sum(tolower(X$source) == "official"),
                  N_s = sum(tolower(X$source) == "survey"),
                  N_i = max(X$i),
                  N_j = length(unique(X$vaccine)),
                  N_t = max(X$t),
                  start_o = min(which(tolower(X$source) == "official")),
                  start_s = min(which(tolower(X$source) == "survey")))

  return(vax_dat)
}

