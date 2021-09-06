
#' Summary methods for 'ic' fit objects
#'
#' Summaries of paramter estimates and MCMC convergence diagnostics (Monte Carlo
#' error, effective sample size, Rhat).
#' @param X Object of type \code{icfit} or \coe{iclist}.
#' @param ... Currently ignored.
#' @param pars Optional character vector specifying a subset of parameters to
#'   display.
#' @param probs Optional numeric vector of probability passed to
#'   \code{\link[stats]{quantile}}.
#' @param digits Number of digits to use for formatting numbers when printing.
#'
#' @return The \code{summary} method returns an object of class
#'   \code{summary.icfit}, which is a matrix of summary statistics and
#'   diagnostics, with attributes storing information for use by the
#'   \code{print} method. The \code{print} method for \code{summary.icfit}
#'   objects is called for its side effect and returns its input.
#'
#' @export
#' @method summary icfit
summary.icfit <- function(X,
                          pars = NULL,
                          probs = c(0.1, 0.5, 0.9),
                          ...,
                          digits = 2){

  if(!is.null(pars)){
    args <- list(object = X$fit, pars = pars, probs = probs)
  } else{
    args <- list(object = X$fit, probs = probs)
  }
  out <- do.call(rstan::summary, args)$summary

  stats <- colnames(out)
  if ("n_eff" %in% stats) {
    out[, "n_eff"] <- round(out[, "n_eff"])
  }

  structure(
    out,
    print.digits = digits,
    class = "summary.icfit"
  )
}


#' @rdname summary.icfit
#' @export
#' @method summary iclist
summary.iclist <- function(X, ...){
  lapply(X, function(dat){ summary(dat, ...) })
}


#' @rdname summary.icfit
#' @export
#' @method print summary.icfit
#'
#' @param X An object of class \code{"summary.icfit"}.
print.summary.icfit <- function(X,
                                digits = max(1, attr(X, "print.digits")),
                                ...){
  cat("\nEstimates:\n")
  sel <- which(colnames(X) %in% c('se_mean', 'n_eff', 'Rhat'))
  has_diagnostic <- length(sel) > 0

  if(has_diagnostic){
    xtemp <- X[, -sel, drop = FALSE]
    colnames(xtemp) <- paste("", colnames(xtemp))
  } else{
    xtemp <- X
  }

  print(format(round(xtemp, digits), nsmall = digits), quote = FALSE, ...)

  if(has_diagnostic){
    cat("\n\nMCMC diagnostics:\n")
    mcse_hat <- format(round(X[, c("se_mean", "Rhat"), drop = FALSE], digits),
                       nsmall = digits)
    n_eff <- format(X[, "n_eff", drop = FALSE], drop0trailing = TRUE)
    print(cbind(mcse_hat, n_eff), quote = FALSE)
    cat("\nFor each parameter, mcse is Monte Carlo standard error,",
        "n_eff is a crude measure of effective sample size,",
        "and Rhat is the potential scale reduction factor on split chains",
        "(at convergence Rhat=1).\n")
  }

  invisible(X)
}


#' @rdname summary.icfit
#' @export
#' @method print summary.iclist
#'
#' @param X An object of class \code{"summary.iclist"}.
print.summary.iclist <- function(X,
                                 digits = max(1, attr(X, "print.digits")),
                                 ...){
  lapply(X, function(dat){ summary(dat, digits, ...) })
}

