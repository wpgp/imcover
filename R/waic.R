
#' Extract  posterior sample of the pointwise log-likelihood from a icfit object.
#'
#' Convenience function to extract a log-likelihood matrix or array from the
#' fitted Stan model in an \code{icfit} object.
#' @param X object of type \code{icfit}
#' @param permuted Logical. Should the draws from each chain be permuted and
#'   merged?. Default is \code{TRUE} to merge.
#' @param ... Not currently used.
#'
#' @return An S x N matrix containing pointwise log-likelihood samples, where S
#'   is the number of samples and N is the number of observations in the data.
#'   If \code{permuted} is \code{FALSE} then an S x N x R array is returned,
#'   where R is the number of chains.
#'
#' @seealso \code{\link[loo]{extract_log_lik}}
#'
#' @aliases log_lik
#' @importFrom rstantools log_lik
#' @importFrom loo extract_log_lik
#' @method log_lik icfit
#' @export
#' @export log_lik
log_lik.icfit <- function(X, permuted = TRUE, ...){

  log_lik <- loo::extract_log_lik(stanfit = X$fit,
                                  parameter_name = 'log_lik',
                                  merge_chains = permuted)

  return(log_lik)
}


#' Widely applicable information criterion (WAIC)
#'
#' Implementation of `waic()` methods from \code{loo} to compute WAIC from the
#' pointwise log-likelihood.
#' @param X Object of type \code{icfit} or \code{iclist}.
#' @param ... Additional arguments. See \code{loo::waic()}.
#' @details For more details see \pkg{loo}. The developers of \code{stan}
#'   recommend LOO-CV using PSIS (as implemented by the \code{loo::loo()}
#'   function) because PSIS provides useful diagnostics as well as effective
#'   sample size and Monte Carlo estimates.
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
#' @seealso \code{\link[loo]{waic}}
#'
#' @aliases waic
#' @importFrom loo waic
#' @method waic icfit
#' @export waic
#' @export
waic.icfit <- function(X, ...){
  LLarray <- log_lik.icfit(X)

  return(loo::waic(LLarray, ...))
}


#' @aliases waic
#' @rdname waic.icfit
#' @importFrom loo waic
#' @method waic iclist
#' @export waic
#' @export
waic.iclist <- function(X, ...){

  lapply(X, FUN = function(dat){ waic(dat, ...) })
}


#' Leave-one-out cross-validation
#'
#' The \code{loo} method for \code{icfit} objects. computes the approximate
#' leave-one-out cross-validation using Pareto smoothed importance sampling
#' (PSIS-LOO CV).
#' @param X Object of type \code{icfit} or \code{iclist}
#' @param ... Additional parameters used by \code{rstan::loo}.
#' @return an object of class \code{loo} with the PSIS results.
#' @seealso \code{\link[loo]{loo}}, \code{\link[rstan]{loo}}
#'
#' @aliases loo
#' @importFrom rstan loo
#' @method loo icfit
#' @export loo
#' @export
loo.icfit <- function(X, ...){

  return(rstan::loo(X$fit, ...))
}


#' @aliases loo
#' @rdname loo.icfit
#' @importFrom rstan loo
#' @method loo iclist
#' @export loo
#' @export
loo.iclist <- function(X, ...){

  lapply(X, FUN = function(dat){ loo(dat, ...) })

}

