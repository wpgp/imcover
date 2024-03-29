% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/posterior_interval.R
\name{ic_coverage}
\alias{ic_coverage}
\alias{ic_coverage.icfit}
\alias{ic_coverage.iclist}
\title{Summarise posterior estimates of immunisation coverage}
\usage{
ic_coverage(
  X,
  object = "posterior",
  stat = c("mean", "median", "sd", "quantile"),
  probs = c(0.025, 0.25, 0.5, 0.75, 0.975)
)

\method{ic_coverage}{icfit}(
  X,
  object = "posterior",
  stat = c("mean", "median", "sd", "quantile"),
  probs = c(0.025, 0.25, 0.5, 0.75, 0.975)
)

\method{ic_coverage}{iclist}(
  X,
  object = "posterior",
  stat = c("mean", "median", "sd", "quantile"),
  probs = c(0.025, 0.25, 0.5, 0.75, 0.975)
)
}
\arguments{
\item{X}{A fitted model object of type \code{icfit} or \code{iclist}.}

\item{object}{Character specifying the name of the object to summarise.
Default is 'posterior'.}

\item{stat}{Character vector of summary statistics to apply.}

\item{probs}{Numeric vector of probabilities with values in \link{0,1} to be used
for \code{stat} "quantile".}
}
\value{
A \code{data.frame} with columns for 'country', 'time', and 'vaccine'
labels and summary statistics, and as many rows as there are 'mu'
observations of immunisation coverage.
}
\description{
Computes the Bayesian posterior summary intervals, sometimes
called 'credible intervals' for the estimated immunisation coverage term
('mu').
}
