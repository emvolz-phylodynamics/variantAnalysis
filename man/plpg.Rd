% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/frequency0.R
\name{plpg}
\alias{plpg}
\title{Performs a simple logistic model fit and make a frequency plot}
\usage{
plpg(
  s1,
  bw = 1,
  cols = c(ancestral = "white", mutant = "black"),
  g = 365/6.5,
  RD = seq(1.1, 1.5, length = 100),
  gentime = 6.5/365,
  timeScaleMethod = c("gentime", "growth")
)
}
\arguments{
\item{s1}{data frame with columns genotype and sample_time}

\item{bw}{bin width for aggregating data}

\item{cols}{Colours for plotting}
}
\description{
genotype has two values: 'ancestral'(or 'wt') and mutant
}
