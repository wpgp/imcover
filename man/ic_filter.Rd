% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter.R
\name{ic_filter}
\alias{ic_filter}
\alias{ic_filter.ic.df}
\alias{ic_filter.icfit}
\alias{ic_filter.iclist}
\title{Filter immunisation coverage records}
\usage{
ic_filter(X, region, country, time, vaccine)

\method{ic_filter}{ic.df}(X, region, country, time, vaccine)

\method{ic_filter}{icfit}(X, country, time, vaccine)

\method{ic_filter}{iclist}(X, region, country, time, vaccine)
}
\arguments{
\item{X}{object of type \code{ic.df}, \code{icfit}, or \code{iclist}}

\item{region}{Character vector of regional groups to select. Not valid for
\code{icfit}.}

\item{country}{Character vector of country identifiers to select}

\item{time}{Vector of time points to select}

\item{vaccine}{Character vector of vaccine records to select, given by
abbreviations (e.g. 'DTP1')}
}
\value{
A new object with type matching \code{X} with the selected records.
}
\description{
Select records from 'ic' objects based on region group, country, time, or
vaccine.
}
\details{
For \code{icfit} and \code{iclist} objects, filtering will only
affect the posterior estimates and prediction (if present). When a
selection parameter is omitted, then all possible values will be selected.
Selection parameters which are not present in \code{X} will be ignored.
}
