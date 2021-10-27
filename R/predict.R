
#' Posterior predictions of vaccination coverage
#'
#' Draw from the posterior predictive distribution.
#' @param X object of type \code{icfit} or \code{iclist}
#' @param country character vector of country codes to predict
#' @param vaccine character vector of vaccine abbreviations to predict
#' @param t Integer of the number of time steps ahead to predict. Default is 2.
#' @param return_ic Logical. Should an \code{icfit} object be returned?
#' @return If \code{return_ic} is \code{FALSE}, then a \code{data.frame} (or a
#'   list of data frames) with posterior samples for predictions, labelled by
#'   'country', 'vaccine', and 'time' columns. Returning an 'ic' object will
#'   modify object \code{X} and add/update with a the 'prediction' element and
#'   return the same type as \code{X}.
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
  n_i <- length(X$labels$lbl_c)
  n_j <- ncol(alpha_j)
  t0 <- ncol(gamma_t)
  max_t <- X$labels$lbl_t[t0]
  nsamples <- nrow(lambda)

  for(tt in (t0 + 1:t)){
    gamma_t <- cbind(gamma_t, matrix(nrow = nrow(gamma_t), ncol = 1))
    gamma_t[, tt] <- rnorm(nsamples, rho_t * gamma_t[, tt-1], sigma_t)
  }

  # modify arrays
  delta_jt <- abind::abind(delta_jt, array(NA, dim = c(nsamples, n_j, t)), along = 3)

  for(ss in 1:nsamples){
    for(tt in (t0 + 1:t)){
      delta_jt[ss, , tt] <- rnorm(n_j, rho_j[ss] * delta_jt[ss, , tt-1], sigma_jt[ss])
    }
  }

  # multi-country models
  if(n_i > 1){
    phi_it <- abind::abind(phi_it, array(NA, dim = c(nsamples, n_i, t)), along = 3)

    for(ss in 1:nsamples){
      for(tt in (t0 + 1:t)){
        phi_it[ss, , tt] <- rnorm(n_i, rho_i[ss] * phi_it[ss, , tt-1], sigma_it[ss])
      }
    }

    # modify array for prediction point
    omega_ijt <- abind::abind(omega_ijt, array(NA, dim = c(nsamples, n_i, n_j, t)), along = 4)

    for(ss in 1:nsamples){
      for(jj in 1:n_j){
        for(tt in (t0 + 1:t)){
          omega_ijt[ss, , jj, tt] <- rnorm(n_i, rho_ij[ss] * omega_ijt[ss, , jj, tt-1], sigma_ijt[ss])
        }
      }
    }
  }


  # combine to generate mu + predictions
  mu_pred <- matrix(NA, nrow = nsamples, ncol = t * n_i * n_j)
  # flat mu_pred index
  mu_id <- expand.grid(ii = 1:n_i, jj = 1:n_j, tt = t0+(1:t))
  mu_id <- mu_id[order(mu_id$ii, mu_id$jj, mu_id$tt), ]

  # generate prediction
  for(ss in 1:nsamples){
    for(idx in 1:nrow(mu_id)){
      ii <- mu_id$ii[idx]
      jj <- mu_id$jj[idx]
      tt <- mu_id$tt[idx]

      if(n_i > 1){
        mu_pred[ss, idx] <- beta_i[ss, ii] + alpha_j[ss, jj] + gamma_t[ss, tt] +
          psi_ij[ss, ii, jj] + phi_it[ss, ii, tt] + delta_jt[ss, jj, tt] +
          omega_ijt[ss, ii, jj, tt]
      } else{
        mu_pred[ss, idx] <- alpha_j[ss, jj] + gamma_t[ss, tt] + delta_jt[ss, jj, tt]
      }
    }
  }

  ## create prediction data frame
  mu <- cbind(mu, mu_pred)
  mu <- t(as.data.frame(mu))
  mu <- invlogit(mu)

  # create an index to the 'mu' full parameter
  mu_id_pred <- mu_id
  # flat mu_pred index
  mu_id <- expand.grid(ii = 1:n_i, jj = 1:n_j, tt = 1:t0)
  mu_id <- mu_id[order(mu_id$ii, mu_id$jj, mu_id$tt), ]
  mu_id <- rbind(mu_id, mu_id_pred)
  names(mu_id) <- c("country", "vaccine", "time")

  # name the indices -- add in the new prediction time points to the labels
  mu_names <- cbind.data.frame(
    'country' = X$labels$lbl_c[mu_id$country],
    'time' = c(X$labels$lbl_t, (max_t + 1:t))[mu_id$time],
    'vaccine' = X$labels$lbl_v[mu_id$vaccine]
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


#' @aliases predict
#' @importFrom stats predict
#' @export predict.iclist
#' @export
predict.iclist <- function(X, country, vaccine, t = 2, return_ic = TRUE){
  out <- lapply(X, FUN = function(fit){
    predict(fit, country, vaccine, t, return_ic)
  })

  if(return_ic){
    class(out) <- list("iclist", class(out))
  }
  return(out)
}



