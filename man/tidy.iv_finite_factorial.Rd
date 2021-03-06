% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iv_finite_factorial.R
\name{tidy.iv_finite_factorial}
\alias{tidy.iv_finite_factorial}
\title{Tidy an iv_finite_factorial object}
\usage{
\method{tidy}{iv_finite_factorial}(x, conf.level = 0.95, ...)
}
\arguments{
\item{x}{An \code{iv_factorial} object produced by a call to
\code{\link[=iv_finite_factorial]{iv_finite_factorial()}}}

\item{conf.level}{The confidence level to use for the confidence
interval. Must be strictly greater than 0 and less than 1.
Defaults to 0.95, which corresponds to a 95 percent confidence
interval.}

\item{...}{Additional arguments. Not used. Needed to match generic
signature only.}
}
\value{
A \code{\link[tibble:tibble]{tibble::tibble()}} with columns:

\item{term}{The name of the effect term.}

\item{estimand}{Which complier effect being
estimated.}

\item{estimate}{The estimated value of the effect.}

\item{ci_1_lower}{Lower bound for the first interval of the
Fieller confidence regresion for the estimate.}

\item{ci_1_upper}{Upper bound for the first interval of the
Fieller confidence regresion for tshe estimate.}

\item{ci_2_lower}{Lower bound for the second interval of the
Fieller confidence regresion for the estimate. Only non-\code{NA}
when the confidence region is disjoint.}
\item{ci_2_upper}{Upper
bound for the second interval of the Fieller confidence
regresion for the estimate. Only non-\code{NA} when the confidence
region is disjoint.}
}
\description{
Tidy summarizes information about the components of a model.
}
\author{
Matt Blackwell
}
