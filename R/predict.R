
#' Posterior predictions of vaccination coverage
#'
#' Draw from the posterior predictive distribution.
#' @param X object of type \code{icfit} or \code{iclist}
#' @param country character vector of country codes to predict
#' @param vaccine character vector of vaccine abbreviations to predict
#' @param t Integer of the number of time steps ahead to predict. Default is 2.
#' @param return_ic Logical. Should an \code{icfit} object be returned?
#' @return If \code{return_ic} is \code{FALSE}, then a \code{data.frame} with
#'   posterior samples for predictions, labelled by 'country', 'vaccine', and
#'   'time' columns. Return and 'ic' object will modify object \code{X} and
#'   add/update a the 'prediction' element.
#'
#' @aliases predict
#' @importFrom stats predict
#' @export predict.icfit
#' @export
predict.icfit <- function(X, country, vaccine, t = 2, return_ic = TRUE){

  # validate inputs
  if(missing(country)){
    country <- X$labels$lbl_c
  } else{
    country <- intersect(country, X$labels$lbl_c)
    stopifnot(length(country) >= 1L)
  }

  if(missing(vaccine)){
    vaccine <- X$labels$lbl_v
  } else{
    vaccine <- intersect(vaccine, X$labels$lbl_v)
    stopifnot(lenth(vaccine) >= 1L)
  }

  stopifnot(t >= 1)
  if(t > 3){
    warning("Predicting beyond 3 time steps is highly uncertain!",
            call. = FALSE, immediate. = TRUE)
  }

  # extract fitted model
  ext <- rstan::extract(X$fit)
  list2env(ext, envir = environment())

  # get key dimensions
  n_i <- ncol(beta_i)
  n_j <- ncol(alpha_j)
  t0 <- ncol(phi_t)
  max_t <- X$labels$lbl_t[t0]
  nsamples <- nrow(beta_i)

  for(tt in (t0 + 1:t)){
    phi_t <- cbind(phi_t, matrix(nrow = nrow(phi_t), ncol = 1))
    phi_t[, tt] <- rnorm(nsamples, rho_t * phi_t[, tt-1], sigma_t)
  }

  # modify arrays
  gamma_it <- abind::abind(gamma_it, array(NA, dim = c(nsamples, n_i, t)), along = 3)
  gamma_jt <- abind::abind(gamma_jt, array(NA, dim = c(nsamples, n_j, t)), along = 3)

  for(ss in 1:nsamples){
    for(tt in (t0 + 1:t)){
      gamma_it[ss, , tt] <- rnorm(n_i, rho_i[ss] * gamma_it[ss, , tt-1], sigma_it)
      gamma_jt[ss, , tt] <- rnorm(n_j, rho_j[ss] * gamma_jt[ss, , tt-1], sigma_jt)
    }
  }

  # modify array for prediction point
  delta_ijt <- abind::abind(delta_ijt, array(NA, dim = c(nsamples, n_i, n_j, t)), along = 4)

  for(ss in 1:nsamples){
    for(jj in 1:n_j){
      for(tt in (t0 + 1:t)){
        delta_ijt[ss, , jj, tt] <- rnorm(n_i, rho_ij[ss] * delta_ijt[ss, , jj, tt-1], sigma_ijt)
      }
    }
  }

  # combine to generate mu + predictions
  mu <- abind::abind(mu, array(NA, dim = c(nsamples, n_i, t, n_j)), along = 3)

  for(ss in 1:nsamples){
    for(ii in 1:n_i){
      for(jj in 1:n_j){
        for(tt in (t0 + 1:t)){
          mu[ss, ii, tt, jj] <- beta_i[ss, ii] + alpha_j[ss, jj] + phi_t[ss, tt] +
            gamma_it[ss, ii, tt] + gamma_jt[ss, jj, tt] + gamma_ij[ss, ii, jj] +
            delta_ijt[ss, ii, jj, tt]
        }
      }
    }
  }

  ## create prediction data frame
  mu <- t(as.data.frame(mu))
  mu <- invlogit(mu)

  # create an index to the 'mu' parameter
  mu_idx <- strsplit(rownames(mu), '.', fixed = T)
  mu_idx <- do.call(rbind.data.frame, mu_idx)
  names(mu_idx) <- c("country", "time", "vaccine")
  mu_idx <- as.data.frame(lapply(mu_idx, function(x) as.numeric(x)))

  # name the indices -- add in the new prediction time points to the labels
  mu_names <- cbind.data.frame(
    'country' = X$labels$lbl_c[mu_idx$country],
    'time' = c(X$labels$lbl_t, (max_t + 1:t))[mu_idx$time],
    'vaccine' = X$labels$lbl_v[mu_idx$vaccine]
  )

  # reverse ratio calculation
  # STILL TO-DO!!!

  # convert to coverage percentage
  mu <- mu * 100

  # apply labels
  mu <- cbind(mu_names, mu)
  rownames(mu) <- NULL

  # prepare data to return
  # subset only prediction times
  mu <- mu[mu$country %in% country &
             mu$time %in% seq(max_t + 1, length.out = t) &
             mu$vaccine %in% vaccine, ]

  if(return_ic){
    X[['prediction']] <- mu
    return(X)

  } else{
    return(mu)
  }
}

