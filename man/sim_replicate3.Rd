% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/power1.R
\name{sim_replicate3}
\alias{sim_replicate3}
\title{Three methods: gam adjusted for imports; simple logistic; hier bayes 
also simulates breakthrough}
\usage{
sim_replicate3(
  MU = lubridate::decimal_date(as.Date("2021-04-18")),
  TFIN,
  RHO0,
  k = 25,
  bs_replicates = 100,
  breakthroughrate = 0.41,
  or_breakthrough = 1.25,
  ncpu = 1,
  ...
)
}
\description{
Three methods: gam adjusted for imports; simple logistic; hier bayes 
also simulates breakthrough
}
