
#' Multi-source immunisation coverage model with Stan
#'
#' @param x Numeric vector of input values.
#' @param y Numberic vector of output values.
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'
#' @name multi_lik_stan
#' @export
multi_lik_stan <- function(x, y, ...) {
  # standata <- list(x = x, y = y, N = length(y))
  out <- rstan::sampling(stanmodels$multi_lik, ...)
  return(out)
}
