
#' Simulate cluster trajectories 
#'
#' SDE model for simple exponential growth, stores cumulative and active infections
#'
#' @param tstart Time of origin
#' @param tfin Simulate until 
#' @param R Expected reproduction number
#' @param Rsd std dev of R 
#' @param non_extinction Minimum cumulative infections required to return trajectory 
#' @return Matrix with simulation output 
sim_cluster_growth <- function(tstart, tfin, R
	, Rsd = .2
	, non_extinction = 0
	, ntries = 10  )
{
	if ( tstart >= tfin ){
		return ( NULL )
	}
	.R <- max( 0, rnorm( 1, R , Rsd ))
	finci = -Inf 
	fn = system.file( 'stochastic_seir0.R', package = 'variantAnalysis' )
	seir_gen <- odin::odin( fn ) 
	seir_sim <- seir_gen( R = .R )
	tries <- 0 
	while( (finci < non_extinction) & (tries < ntries) ) {
		X = seir_sim$run( seq( tstart, tfin, by = 1) ) 
		finci <- tail( X[, 'cI'], 1 )
		print( finci )
		tries <- tries + 1
	}
	
	cbind( X, R = .R )
}

.sim_cluster_growth <- function(tstart, tfin, R
	, Rsd = .2
	, non_extinction = 0
	, ntries = 10  )
{
	#~ 	tstart = 1 
	#~ 	tfin = 30 
	#~ 	R = 1.3 
	if ( tstart >= tfin ){
		return ( NULL )
	}
	.R <- max( 0, rnorm( 1, R , Rsd ))
	finci = -Inf 
	fn = system.file( 'stochastic_seir1.R', package = 'variantAnalysis' )
	seir_gen <- odin::odin( fn ) 
	seir_sim <- seir_gen( R = .R , tini = tstart, Requil = R )
	tries <- 0 
	while( (finci < non_extinction) & (tries < ntries) ) {
		X = seir_sim$run( seq( tstart, tfin, by = 1) ) 
		finci <- tail( X[, 'cI'], 1 )
		print( finci )
		ntries <- ntries + 1
	}
	
	cbind( X, R = .R )
}

#~ sim_cluster_growth( 0, 30 , 1.2, non_extinction = 5 / .1 / .5  )


#' Simulate time and number of importation events for clusters of two lineages
#' 
#' @param E_imports0 Expected number of imports for variant 0 which initiate cluster
#' @param mu Mean time of import
#' @param sigma std dev of importation time
#' @return $variant and $ancestral contain times of import 
sim_importation_time <- function(E_imports0, E_imports1, mu, sigma ) {
	n0 <- rpois( 1, E_imports0 )
	n1 <- rpois( 1, E_imports1 )
	list(
		variant = rnorm( n0, mu, sigma )
		, ancestral = rnorm( n1, mu, sigma )
	)
}


#' Simulate sampling times given cluster trajectory and sample rate
#'
#' Note sampling rate decays to zero from tfin - 12 days to tfin
#'
#' @param X output of sim_cluster_growth
#' @param rho0 probability that case is sampled 
#' @param rho1 probability that case is diagnosed 
#' @return vector of sample times 
sim_sampling_cluster <- function( X,  rho0 = .10 , rho1 = .25 ) 
{
	N <- floor( tail(X[, 'cI'], 1) )
	n <- rbinom( 1, N, rho0*rho1 )
	w <- pmax(0, as.vector( X[, 'I'] ) )
	lb <- min( 11, nrow(X)-1)
	taili <- (length(w ) - lb):length(w)
	w[taili] <- w[taili] *  exp( -(0:lb) / 3 )  # seq( 1-1/12, 0, by = -1/12 )
	if ( sum( w ) > 0 ){
		w  <- w / sum( w ) 
	} else{
		w <- rep(1 , length( w) )
	}
	if ( any( is.na( w ))) {
		browser()  
		return( NULL )
	}
	ti <- sample.int( nrow(X)
	 , size = n 
	 , prob = w
	 , replace = TRUE
	)
	as.vector( X[ti, 1] )
}


#' Simulate entire pipeline: importation, growth, sampling 
sim_replicate <- function( 
 E_imports0 = 650 # cf coguk/g2-adf3.rds
 , E_imports1 = 650
 , mu = lubridate::decimal_date( as.Date( "2021-04-18")) #cf coguk/g2-adf3.rds, based on B117
 , sigma = 12.9 / 365 # std dev of cluster origin (yrs), g2-adf3.rds
 , tfin = lubridate::decimal_date( as.Date( "2021-04-18")) + 1/12
 , Rvariant = 1.6*1. #cf coguk/b.1.617/g3.1
 , Rancestral = 1.
 , Rsd = .2 # std dev of initial R in clusters, cf coguk/b.1.617/g3.1
 , rho0 = 0.9 # proportion sequenced  
 , rho1 = 0.5 # proportion diagnosed 
)
{ 
	library( lubridate )
	rho <- rho0 * rho1 
	nonext <- 5 / rho 
	
	itimes <- sim_importation_time(E_imports0, E_imports1, mu, sigma )
	
	## days
	.tfin <- as.numeric( as.Date( date_decimal( tfin ) ) - as.Date( '2021-01-01' ) )
	
	ctraj0 <- lapply( itimes$variant, function(tt){
		.tstart = as.numeric( as.Date( date_decimal( tt ) ) - as.Date( '2021-01-01' ) )
		sim_cluster_growth(.tstart, .tfin, Rvariant, Rsd= Rsd,  non_extinction = 10 )
	})
	samp0 <- do.call( rbind, lapply( seq_along(ctraj0), function(i){
		if ( is.null( ctraj0[[i]] ) )
			return( NULL )
		ssc =  sim_sampling_cluster( ctraj0[[i]],  rho0 = rho0, rho1 = rho1 ) 
		if ( length( ssc ) == 0 )
			return( NULL )
		data.frame( cluster = paste(sep='.', 'variant', i )
		 , sample_time = ssc
		 , lineage = 'variant'
		 , stringsAsFactors = FALSE
		)
	}))
	
	ctraj1 <- lapply( itimes$ancestral, function(tt){
		.tstart = as.numeric( as.Date( date_decimal( tt ) ) - as.Date( '2021-01-01' ) )
		sim_cluster_growth(.tstart, .tfin, Rancestral, Rsd= Rsd, non_extinction = 10 )
	})
	samp1 <- do.call( rbind, lapply( seq_along(ctraj1), function(i){
		if ( is.null( ctraj1[[i]] ) )
			return( NULL )
		ssc =  sim_sampling_cluster( ctraj1[[i]],  rho0 = rho0, rho1 = rho1 ) 
		if ( length( ssc ) == 0 )
			return( NULL )
		data.frame( cluster = paste(sep='.', 'ancestral', i )
		 , sample_time = ssc
		 , lineage = 'ancestral'
		 , stringsAsFactors = FALSE
		)
	}))
	
	rbind( samp0, samp1 )	
}

#~ library( lubridate )
#~ o = sim_replicate( 
#~  E_imports0 = 650 # cf coguk/g2-adf3.rds
#~  , E_imports1 = 650
#~  , mu = decimal_date( as.Date( "2021-04-18")) #cf coguk/g2-adf3.rds, based on B117
#~  , sigma = 12.9 / 365 # std dev of cluster origin (yrs), g2-adf3.rds
#~  , tfin = decimal_date( as.Date( "2021-04-18")) + 1/12
#~  , Rvariant = 1.6*0.9 #cf coguk/b.1.617/g3.1
#~  , Rancestral = .9
#~  , Rsd = .2 # std dev of initial R in clusters, cf coguk/b.1.617/g3.1
#~  , rho0 = 0.9 # proportion sequenced  
#~  , rho1 = 0.5 # proportion diagnosed TODO TODO 
#~ )


#~ #' Inference for simulated samples

#~ library( mlogit ) 
#~ minsize <- 20 
#~ tc <- table( o$cluster )
#~ keepc <- names(tc) [ tc >= minsize ] 
#~ o1 <- o[ o$cluster %in% keepc , ]

#~ m = glm( (lineage=='variant') ~ sample_time , family = binomial( link = 'logit' ) , data = o ); summary( m )
#~ m = mgcv::gam( as.factor( cluster )  ~ sample_time * lineage , data = o1 , family = mgcv::multinom(K = length(unique( o1$cluster)))  )

