
#' Widely applicable information criterion (WAIC)
#'
#' Implementation of `waic()` methods from \code{loo} to compute WAIC from the
#' pointwise log-likelihood.
#' @param X Object of type \code{icfit} or \code{iclist}.
#' @param pars Character names of parameters to extract. Default is the
#'   log-likelihood.
#' @param ... Additional arguments. See \code{loo::loo()} and
#'   \code{loo::waic()}.
#' @param save_psis Should the "psis" object created internally by loo() be
#'   saved in the returned object? Default is \code{FALSE}.
#' @param cores The number of cores to use for parallelization. This defaults to
#'   the option mc.cores which can be set for an entire R session by
#'   options(mc.cores = NUMBER).
#' @param ... Additional arguments passed to \code{loo::loo()}.
#' @details The developers of \code{stan} recommend LOO-CV using PSIS (as
#'   implemented by the \code{loo::loo()} function) because PSIS provides useful
#'   diagnostics as well as effective sample size and Monte Carlo estimates.
#'
#' @return A named list (of class `c("waic", "loo")`) with components:
#'
#' \describe{
#'  \item{`estimates`}{
#'  A matrix with two columns (`"Estimate"`, `"SE"`) and three
#'  rows (`"elpd_waic"`, `"p_waic"`, `"waic"`). This contains
#'  point estimates and standard errors of the expected log pointwise predictive
#'  density (`elpd_waic`), the effective number of parameters
#'  (`p_waic`) and the information criterion `waic` (which is just
#'  `-2 * elpd_waic`, i.e., converted to deviance scale).
#'  }
#'  \item{`pointwise`}{
#'  A matrix with three columns (and number of rows equal to the number of
#'  observations) containing the pointwise contributions of each of the above
#'  measures (`elpd_waic`, `p_waic`, `waic`).
#'  }
#' }
#'
#' @references
#' Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation and
#' widely application information criterion in singular learning theory.
#' *Journal of Machine Learning Research* **11**, 3571-3594.
#'
#' @seealso \code{[loo](waic)}
#'
#' @aliases waic
#' @importFrom loo waic waic.function waic.matrix is.waic
#' @export
#' @method waic icfit
waic.icfit <- function(X, ...){
  LLarray <- loo::extract_log_lik(stanfit = X$fit,
                                  parameter_name = 'log_lik',
                                  merge_chains = TRUE)
  loo::waic(LLarray, ...)
}


#' @aliases waic
#' @importFrom loo waic waic.function waic.matrix is.waic
#' @export
#' @method waic iclist
waic.iclist <- function(X, ...){

  lapply(X, FUN = function(dat){ waic(dat, ...) })
}


#' @aliases loo
#' @importFrom loo loo loo.function loo.matrix is.loo
#' @export
#' @method loo icfit
#' @rdname waic.icfit
loo.icfit <- function(X,
                      pars = "log_lik",
                      ...,
                      save_psis = FALSE,
                      cores = getOption("mc.cores", 1)){

  stopifnot(length(pars) == 1L)

  LLarray <- loo::extract_log_lik(stanfit = X$fit,
                                  parameter_name = pars,
                                  merge_chains = FALSE)

  r_eff <- loo::relative_eff(x = exp(LLarray), cores = cores)

  loo::loo.array(LLarray,
                 r_eff = r_eff,
                 cores = cores,
                 save_psis = save_psis)
}


#' @aliases loo
#' @importFrom loo loo loo.function loo.matrix is.loo
#' @export
#' @method loo iclist
#' @rdname waic.icfit
loo.iclist <- function(X,
                       pars = "log_lik",
                       ...,
                       save_psis = FALSE,
                       cores = getOption("mc.cores", 1)){

  lapply(X, FUN = function(dat){ loo(dat, pars, ..., save_psis, cores) })
}

