% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/power.R
\name{sim_powerDetectionPrevalence_logisticGrowth}
\alias{sim_powerDetectionPrevalence_logisticGrowth}
\title{Simple power calculation for detecting growth of variant with constant selection coefficient.}
\usage{
sim_powerDetectionPrevalence_logisticGrowth(
  selcoef = 0.25,
  nperday = 10,
  alpha = 0.05,
  Tg = 6.5,
  mindays = 31,
  maxdays = 240,
  init_freq = 1e-04,
  nreplicates = 5
)
}
\arguments{
\item{selcoef}{selection coefficient (per generation)}

\item{nperday}{samples per day within region}

\item{alpha}{Signif threshold for detection}

\item{Tg}{generation time (days)}

\item{mindays}{min days to search for possible detection}

\item{maxdays}{max days to search for possible detection}

\item{init_freq}{Initial frequency of variant (impacts computation, will not usually need to change this)}

\item{nreplicates}{Results averaged over this many replicates}
}
\description{
This does not account for possible confounders such as importation, changing Rt, outbreaks etc
}
