
#' Multi-source immunisation coverage model with Stan
#'
#' @param X Object of \code{ic.df} for analysis
#' @param prior_lambda Scale parameter for normal prior on source-specific
#'   intercepts. See details.
#' @param prior_sigma Scale parameter for the truncated Cauchy prior on
#'   source-specific standard deviations. See details.
#' @param upper_sigma Numeric value (or vector) setting the upper bounds on the
#'   source-specific scale parameter. See details.
#' @param lower_sigma Numeric value (or vector) setting the lower bounds on the
#'   source-specific scale parameter. See details.
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
#' @details The default priors for the source-specific intercepts (lambda) and
#'   standard deviations (sigma) are given by normal and truncated Cauchy
#'   distributions, respectively. Future version may allow users to specify
#'   different distributions. Currently, users may only specify the scale
#'   parameter of these distribution. By default, lambda for all sources is
#'   given a value of 0.5, i.e. a prior distribution of \code{normal(0, 0.5)}.
#'
#'   For sigma, the default prior varies by source. While administrative and
#'   official data are set to a value of 2, survey data is restricted to 0.2,
#'   reflecting a prior belief in more accurate measurements. These parameters
#'   are used in a half-Cauchy distribution (e.g. \code{cauchy(0, 0.2)}).
#'
#'   In addition, the upper and lower bounds of the sigma parameter may be set
#'   to define the range of the possible values. In general, the lower bounds
#'   should always be zero and cannot be negative. Similar to the differences in
#'   the prior distributions, an upper bound of 0.4 is placed on 'survey'
#'   estimates. In the absence of a user-defined upper-bound, these are set to
#'   100.
#'
#'   Users setting \code{prior_lambda}, \code{prior_sigma}, or
#'   \code{upper_sigma} should note that the length of values specified must be
#'   1 or match the number of unique sources in \code{X}. When only one value is
#'   present, it will be applied to all sources. When a vector is provided, the
#'   order of values must match the order of sources as given by
#'   \code{list_sources(X)}.
#'
#' @name ic_fit
#' @export
ic_fit <- function(X,
                   prior_lambda = NULL, prior_sigma = NULL,
                   upper_sigma = NULL, lower_sigma = NULL,
                   region = TRUE, verbose = TRUE, ...){
  UseMethod("ic_fit")
}


#' @name ic_fit
#' @export
ic_fit.ic.df <- function(X,
                         prior_lambda = NULL, prior_sigma = NULL,
                         upper_sigma = NULL, lower_sigma = NULL,
                         region = TRUE, verbose = TRUE, ...){

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
    out <- multi_lik_stan(X, prior_lambda, prior_sigma,
                          upper_sigma, lower_sigma, verbose, ...)
  } else{
    # main processing loop
    out <- lapply(regions, function(r){
      if(verbose) print(r)
      # subset
      dat <- X[X[[get_attr(X, 'region')]] == r, ]
      # calculate
      fit <- multi_lik_stan(dat, prior_lambda, prior_sigma,
                            upper_sigma, lower_sigma, verbose, ...)

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
#' @param prior_lambda Scale parameter for half-normal prior on source-specific
#'   intercepts. See details.
#' @param prior_sigma Scale parameter for half-Cauchy prior on source-specific
#'   standard deviations. See details.
#' @param upper_sigma Numeric value (or vector) setting the upper bounds on the
#'   source-specific scale parameter. See details.
#' @param lower_sigma Numeric value (or vector) setting the lower bounds on the
#'   source-specific scale parameter. See details.
#' @param verbose Logical. Should messages be displayed? Default is \code{TRUE}.
#' @param ... Arguments passed to \code{rstan::sampling} (e.g. iter, chains).
#' @return An object of class \code{stanfit} returned by \code{rstan::sampling}.
#'
#' @name multi_lik_stan
#' @keywords internal
multi_lik_stan <- function(X,
                           prior_lambda = NULL, prior_sigma = NULL,
                           upper_sigma = NULL, lower_sigma = NULL,
                           verbose = TRUE, ...) {
  if(!is.ic_data(X)) stop("Please provide valid 'ic' data.")

  # data prep
  vax_data <- ic_to_stan(X)
  lbls <- vax_data[[2]]
  vax_data <- vax_data[[1]]

  # flat mu index
  mu_id <- expand.grid(ii = 1:vax_data$N_i, jj = 1:vax_data$N_j, tt = 1:vax_data$N_t)
  mu_id <- mu_id[order(mu_id$ii, mu_id$jj, mu_id$tt), ]
  # lookup - what record in 'mu' for observed coverage 'y'
  mu_lookup <- match(paste(vax_data$i, vax_data$j, vax_data$t), paste(mu_id$ii, mu_id$jj, mu_id$tt))

  # check/set priors
  if(!is.null(prior_lambda)){
    if(length(prior_lambda) == 1){
      vax_data$prior_lambda <- rep(prior_lambda, vax_data$nsources)
    } else{
      if(length(prior_lambda) != vax_data$nsources){
        stop("Priors for lambda must match the number of data sources", call. = FALSE)
      }
      vax_data$prior_lambda <- prior_lambda
    }
  } else{  # defaults
    vax_data$prior_lambda <- rep(0.5, vax_data$nsources)
  }

  if(!is.null(prior_sigma)){
    if(length(prior_sigma) == 1){
      vax_data$prior_sigma <- rep(prior_sigma, vax_data$nsources)
    } else{
      if(length(prior_sigma) != vax_data$nsources){
        stop("Priors for sigma must match the number of data sources", call. = FALSE)
      }
      vax_data$prior_sigma <- prior_sigma
    }
  } else{  # defaults
    prior_sigma <- rep(2, vax_data$nsources)
    prior_sigma[grepl('survey', list_sources(X), fixed = T)] <- 0.2
    vax_data$prior_sigma <- prior_sigma
  }

  if(!is.null(upper_sigma)){
    if(length(upper_sigma) == 1){
      vax_data$U_sigma <- rep(upper_sigma, vax_data$nsources)
    } else{
      if(length(upper_sigma) != vax_data$nsources){
        stop("Priors for sigma must match the number of data sources", call. = FALSE)
      }
      vax_data$U_sigma <- upper_sigma
    }
  } else{  # defaults
    upper_sigma <- rep(100, vax_data$nsources)
    upper_sigma[grepl('survey', list_sources(X), fixed = T)] <- 0.4
    vax_data$U_sigma <- upper_sigma
  }

  if(!is.null(lower_sigma)){
    if(length(lower_sigma) == 1){
      vax_data$L_sigma <- rep(lower_sigma, vax_data$nsources)
    } else{
      if(length(lower_sigma) != vax_data$nsources){
        stop("Priors for sigma must match the number of data sources", call. = FALSE)
      }
      vax_data$L_sigma <- lower_sigma
    }
  } else{
    vax_data$L_sigma <- rep(0, vax_data$nsources)
  }

  if(any(vax_data$L_sigma >= vax_data$U_sigma)){
    stop("Invalid bounds on sigma", call. = FALSE)
  }

  # update model data with mu indexing
  vax_data$ii <- mu_id$ii
  vax_data$jj <- mu_id$jj
  vax_data$tt <- mu_id$tt
  vax_data$mu_lookup <- mu_lookup

  # call stan model
  out <- rstan::sampling(stanmodels$multi_lik_v2,
                         data = vax_data,
                         show_messages = verbose,
                         ...)

  stopifnot(out@mode == 0) # check for model fitting

  # calculate prediction ('mu')
  posterior <- t(as.data.frame(out, 'mu'))
  posterior <- invlogit(posterior)

  # name the indices
  mu_names <- cbind.data.frame(
    'country' = lbls$lbl_c[mu_id$ii],
    'time' = lbls$lbl_t[mu_id$tt],
    'vaccine' = lbls$lbl_v[mu_id$jj]
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
              'model' = 'multi-likelihood',
              'posterior' = posterior,
              'data' = X,  # vax_data[[1]],
              'labels' = lbls,
              'numerator' = attr(X, 'numerator'),
              'denominator' = attr(X, 'denominator'))

  class(out) <- list('icfit', class(out))

  return(out)
}



#' Multi-source immunisation coverage model with Stan
#'
#' @description Sampling from a single likelihood model using a random effect
#'   for data source.
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
#' @details The current version of \code{imcover} does not allow for users to
#'   modify the prior distributions in the single likelihood model. Future
#'   versions of the package may enable this functionality. Lambda is assigned a
#'   prior distribution of \code{normal(0, 1)} and the sigma parameters are
#'   assigned half-Cauchy priors with a scale of 2.
#'
#' @name ic_fit_single
#' @export
ic_fit_single <- function(X, region = TRUE, verbose = TRUE, ...){
  UseMethod("ic_fit_single")
}


#' @name ic_fit
#' @export
ic_fit_single.ic.df <- function(X, region = TRUE, verbose = TRUE, ...){

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
    out <- single_lik_stan(X, verbose, ...)
  } else{
    # main processing loop
    out <- lapply(regions, function(r){
      if(verbose) print(r)
      # subset
      dat <- X[X[[get_attr(X, 'region')]] == r, ]
      # calculate
      fit <- single_lik_stan(dat, verbose, ...)

      return(fit)
    })
    names(out) <- regions
    class(out) <- list("iclist", class(out))
  }

  return(out)
}


#' Multi-source immunisation coverage model with with Stan
#'
#' @description Sampling from a single likelihood model using a random effect
#'   for data source.
#' @param X Object of \code{ic.df} for analysis
#' @param verbose Logical. Should messages be displayed? Default is \code{TRUE}.
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'
#' @name single_lik_stan
#' @keywords internal
single_lik_stan <- function(X, verbose = TRUE, ...) {
  if(!is.ic_data(X)) stop("Please provide valid 'ic' data.")

  # data prep
  vax_data <- ic_to_stan(X)
  lbls <- vax_data[[2]]
  vax_data <- vax_data[[1]]

  # flat mu index
  mu_id <- expand.grid(ii = 1:vax_data$N_i, jj = 1:vax_data$N_j, tt = 1:vax_data$N_t)
  mu_id <- mu_id[order(mu_id$ii, mu_id$jj, mu_id$tt), ]
  # lookup - what record in 'mu' for observed coverage 'y'
  mu_lookup <- match(paste(vax_data$i, vax_data$j, vax_data$t), paste(mu_id$ii, mu_id$jj, mu_id$tt))

  # # check/set priors
  # if(!is.null(prior_lambda)){
  #   if(length(prior_lambda) == 1){
  #     vax_data$prior_lambda <- rep(prior_lambda, vax_data$nsources)
  #   } else{
  #     if(length(prior_lambda) != vax_data$nsources){
  #       stop("Priors for lambda must match the number of data sources", call. = FALSE)
  #     }
  #     vax_data$prior_lambda <- prior_lambda
  #   }
  # } else{  # defaults
  #   vax_data$prior_lambda <- rep(0.5, vax_data$nsources)
  # }
  #
  # if(!is.null(prior_sigma)){
  #   if(length(prior_sigma) == 1){
  #     vax_data$prior_sigma <- rep(prior_sigma, vax_data$nsources)
  #   } else{
  #     if(length(prior_sigma) != vax_data$nsources){
  #       stop("Priors for sigma must match the number of data sources", call. = FALSE)
  #     }
  #     vax_data$prior_sigma <- prior_sigma
  #   }
  # } else{  # defaults
  #   prior_sigma <- rep(2, vax_data$nsources)
  #   prior_sigma[grepl('survey', list_sources(X), fixed = T)] <- 0.2
  #   vax_data$prior_sigma <- prior_sigma
  # }
  #
  # if(!is.null(upper_sigma)){
  #   if(length(upper_sigma) == 1){
  #     vax_data$U_sigma <- rep(upper_sigma, vax_data$nsources)
  #   } else{
  #     if(length(upper_sigma) != vax_data$nsources){
  #       stop("Priors for sigma must match the number of data sources", call. = FALSE)
  #     }
  #     vax_data$U_sigma <- upper_sigma
  #   }
  # } else{  # defaults
  #   upper_sigma <- rep(100, vax_data$nsources)
  #   upper_sigma[grepl('survey', list_sources(X), fixed = T)] <- 0.4
  #   vax_data$U_sigma <- upper_sigma
  # }
  #
  # if(!is.null(lower_sigma)){
  #   if(length(lower_sigma) == 1){
  #     vax_data$L_sigma <- rep(lower_sigma, vax_data$nsources)
  #   } else{
  #     if(length(lower_sigma) != vax_data$nsources){
  #       stop("Priors for sigma must match the number of data sources", call. = FALSE)
  #     }
  #     vax_data$L_sigma <- lower_sigma
  #   }
  # } else{
  #   vax_data$L_sigma <- rep(0, vax_data$nsources)
  # }
  #
  # if(any(vax_data$L_sigma >= vax_data$U_sigma)){
  #   stop("Invalid bounds on sigma", call. = FALSE)
  # }

  # update model data with mu indexing
  vax_data$ii <- mu_id$ii
  vax_data$jj <- mu_id$jj
  vax_data$tt <- mu_id$tt
  vax_data$mu_lookup <- mu_lookup

  # call stan model
  out <- rstan::sampling(stanmodels$single_lik,
                         data = vax_data,
                         show_messages = verbose,
                         ...)

  stopifnot(out@mode == 0) # check for model fitting

  # calculate prediction ('mu')
  posterior <- t(as.data.frame(out, 'mu'))
  posterior <- invlogit(posterior)

  # name the indices
  mu_names <- cbind.data.frame(
    'country' = lbls$lbl_c[mu_id$ii],
    'time' = lbls$lbl_t[mu_id$tt],
    'vaccine' = lbls$lbl_v[mu_id$jj]
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
              'model' = 'single-likelihood',
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
                  nsources = max(X$s),
                  source = X$s,
                  sizes = as.numeric(table(X$source)),
                  N_i = max(X$i),
                  N_j = length(unique(X$vaccine)),
                  N_t = max(X$t)
                 )

  return(list(vax_dat, lbls))
}

