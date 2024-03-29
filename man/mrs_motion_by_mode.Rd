% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_motion_by_mode.R
\name{mrs_motion_by_mode}
\alias{mrs_motion_by_mode}
\alias{mrs_motion_dmsfnj}
\alias{mrs_motion_dhnj}
\alias{mrs_motion_nhnj}
\alias{mrs_motion_rwsipnj}
\title{Calculate mode-by-site motion-response matrices}
\usage{
mrs_motion_dmsfnj(mutants)

mrs_motion_dhnj(mutants)

mrs_motion_nhnj(mutants)

mrs_motion_rwsipnj(mutants)
}
\arguments{
\item{mutants}{A tibble of single-point mutants generated using `generate_mutants`}
}
\value{
a response matrix of the form \eqn{R_{nj}} (response mode is \code{n}, mutated site is \code{j})
}
\description{
Calculate mode-by-site motion-response matrices
}
\details{
`mrs_motion_dmsfnj()` calculates the change of msf along mode n averaged over mutations at site j.

`mrs_motion_dhnj()` calculates change of entropy contribution of mode n averaged over mutations at site j.

`mrs_motion_nhnj()` calculates conservation score nh for mode n averaged over mutations at site j.

`mrs_motion_rwsipnj()` calculates rwsip for mode n averaged over mutations at j
}
