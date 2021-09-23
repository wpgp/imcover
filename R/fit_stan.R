
#' Multi-source immunisation coverage model with Stan
#'
#' @param X Object of \code{ic.df} for analysis
#' @param region Logical. Should region-specific models be generated? Default
#'   is \code{TRUE}.
#' @param verbose Logical. Should messages be displayed? Default is \code{TRUE}.
#' @param ... Arguments passed to \code{rstan::sampling} (e.g. iter, chains).
#' @return An object of class \code{icfit} for single region models, or
#'   \code{iclist} for multiple regions. The \code{icfit} object contains a
#'   \code{stanfit} object returned by \code{rstan::sampling} with the fitted
#'   object along with the posterior samples, data used in the model, and other
#'   attributes related to the fitting. An \code{iclist} object is an extension
#'   of a list containing multiple \code{icfit} objects.
#'
#' @details
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
  vax_data <- ic_to_stan(X)
  lbls <- vax_data[[2]]

  # call stan model
  out <- rstan::sampling(stanmodels$multi_lik,
                         data = vax_data[[1]],
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
    'country' = lbls$lbl_c[mu_idx$country],
    'time' = lbls$lbl_t[mu_idx$time],
    'vaccine' = lbls$lbl_v[mu_idx$vaccine]
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

  out <- list('fit' = out,
              'posterior' = posterior,
              'data' = X,  # vax_data[[1]],
              'labels' = lbls,
              'numerator' = attr(X, 'numerator'),
              'denominator' = attr(X, 'denominator'))

  class(out) <- list('icfit', class(out))

  return(out)
}


#' Prepare 'ic' data for Stan model
#'
#' @param X Object of \code{ic.df} to be converted
#' @return A list with named objects suitable for use as Stan data.
#'   Default is \code{TRUE}.
#' @keywords internal
ic_to_stan <- function(X){
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
  f_v <- factor(X$vaccine)
  f_c <- factor(X$country)
  f_t <- factor(X$time)
  f_s <- factor(X$source)

  X$v <- as.numeric(f_v)
  X$i <- as.numeric(f_c)
  X$t <- as.numeric(f_t) # temporal trend per country
  X$s <- as.numeric(f_s)

  # labels for unique values
  lbl_v <- levels(f_v)
  lbl_c <- levels(f_c)
  lbl_t <- as.numeric(levels(f_t))
  lbl_s <- levels(f_s)

  lbls <- list('lbl_v' = lbl_v, 'lbl_c' = lbl_c,
               'lbl_t' = lbl_t, 'lbl_s' = lbl_s)

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

  return(list(vax_dat, lbls))
}

