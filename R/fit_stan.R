
#' Multi-source immunisation coverage model with Stan
#'
#' @param X Object of \code{ic.df} for analysis
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'
#' @name multi_lik_stan
#' @export
multi_lik_stan <- function(X, ...) {
  if(!is.ic_data(X)) stop("Please provide valid 'ic' data.")

  # check if a ratio adjust was performed
  ratioed <- attr(X, 'ratio')
  # check data
  X <- X[!is.na(X[[get_attr(X, 'coverage')]]), ]
  X <- X[!is.na(X[[get_attr(X, 'region')]]), ]

  # get regions
  regions <- unique(X[[get_attr(X, 'region')]])
  regions <- regions[!is.na(regions)]

  # main processing loop
  out <- lapply(regions, function(r){
    # subset
    dat <- X[X[[get_attr(X, 'region')]] == r, ]

    # data prep
    vax_data <- ic_to_stan(dat, make_index = TRUE)

    fit <- rstan::sampling(stanmodels$multi_lik,
                           data = vax_data
                           ...)
    return(fit)
  })
  names(out) <- regions

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
  nn <- get_attr(X, ic_core())
  X <- data.frame(X)
  idx <- match(nn, names(X))
  names(X)[na.omit(idx)] <- names(nn)[which(!is.na(idx))]

  # sort
  X <- X[order(X$source, X$country, X$vaccine, X$time), ]

  # modify outcome
  X$coverage <- X$coverage / 100
  X$cov.logit <-  log(X$coverage / (1 - X$coverage))

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
