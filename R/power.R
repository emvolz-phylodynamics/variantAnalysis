#' Simple power calculation for detecting growth of variant with constant selection coefficient. 
#'
#' This does not account for possible confounders such as importation, changing Rt, outbreaks etc 
#'
#' @param selcoef selection coefficient (per generation)
#' @param nperday samples per day within region
#' @param alpha Signif threshold for detection
#' @param Tg generation time (days)
#' @param mindays min days to search for possible detection
#' @param maxdays max days to search for possible detection
#' @param init_freq Initial frequency of variant (impacts computation, will not usually need to change this)
#' @param nreplicates Results averaged over this many replicates 
#' @export
sim_powerDetectionPrevalence_logisticGrowth <- function(selcoef = 0.25
 , nperday = 10 
 , alpha = 0.05 
 , Tg = 6.5 
 , mindays = 31
 , maxdays = 240 
 , init_freq = 1e-4 
 , nreplicates = 5 
)
{
	s <- selcoef / Tg 
	
	t0 = log( init_freq / .5 ) / s
	taxis <- (t0):(-t0)
	taxis <- taxis[1:maxdays]
	
	.f <- function(day){
		exp( day * s  ) / ( 1 + exp(  day * s ) )
	}
	f <- .f( taxis )
	
	testindices <- rep(  mindays:maxdays , each = nreplicates )
	pvals = sapply( testindices , function(ub){
		sampf <- rep( f[1:ub], each = nperday )
		.nperday = rpois(ub, lambda = nperday)
		#sampf <- rep( f[1:ub], each = .nperday  )
		sampf = do.call( c, lapply( seq_along( .nperday ), function(i) rep( f[i], .nperday[i] ) ) )
		#sampdays <- rep( taxis[1:ub], each = nperday )
		sampdays <- do.call( c, lapply( seq_along( .nperday ), function(i) rep( taxis[i], .nperday[i] ) ) )
		y = rbinom( length( sampf ) , size = 1, prob = sampf )
		m = glm( y ~ day , data= data.frame( day = sampdays, y = y ), family = binomial( link = 'logit' ) )
		pval = summary( m )$coefficients[2,4 ]
		pval
	})
	logitpvals <- log( pvals) - log( 1 - pvals )
	
	logitpvals1 <- logitpvals 
	logitpvals1[ is.infinite( logitpvals ) ] <- NA 
	#ss = splines2::iSpline( logitpvals1 )
	testtaxis = taxis[ testindices ]
	m = mgcv::gam( logitpvals1 ~ s( testtaxis, bs='cr' ), family= gaussian()  ) 
	of = function( day )
		(predict( m,  newdata = data.frame ( testtaxis=day )) - qnorm( alpha/2 ))^2
	daystar = optimise( of , interval = range( testtaxis ) )$minimum 
	
	rv = .f( daystar )
	print( rv )
	rv
}
