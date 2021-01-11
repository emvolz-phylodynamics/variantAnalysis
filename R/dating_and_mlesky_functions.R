datetree <- function(mltr, ncpu = 6)
{
  sts <- setNames( civmd$sample_time[  match( mltr$tip.label , civmd$central_sample_id ) ], mltr$tip.label )
  tr <- di2multi(mltr, tol = 1e-05)
  tr = unroot(multi2di(tr))
  tr$edge.length <- pmax(1/29000/5, tr$edge.length)
  dater(unroot(tr), sts[tr$tip.label], s = 29000, omega0 = meanrate, numStartConditions = 0, meanRateLimits = c(meanrate, meanrate + 1e-6), ncpu = ncpu)
}



bootrep <- function(mltr, ...)
{
  td = datetree( mltr , ... ) 
  res = diff( range( epiweek( date_decimal( td$sts )  )  )  ) + 1
  res <- res * 2
  sg = mlskygrid( td, tau = NULL, tau_lower=.001, tau_upper = 10 , sampleTimes = td$sts , res = 12, ncpu = 6)
  res = with( sg, approx( time, ne, rule = 1, xout = taxis )$y )
  print( res ) 
  res 
}

