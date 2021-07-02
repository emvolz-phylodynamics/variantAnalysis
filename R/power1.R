
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
	, Rsd = .15
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
	fn = system.file( 'stochastic_seir0.R', package = 'variantAnalysis' )
	seir_gen <- odin::odin( fn ) 
	seir_sim <- seir_gen( R = .R )
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
	taili <- (length(w ) - 11):length(w)
	w[taili] <- w[taili] *  exp( -(0:11) / 3 )  # seq( 1-1/12, 0, by = -1/12 )
	w  <- w / sum( w )  
	ti <- sample.int( nrow(X)
	 , size = n 
	 , prob = w
	 , replace = TRUE
	)
	as.vector( X[ti, 1] )
}


#' Simulate entire pipeline: importation, growth, sampling 
sim_replicate <- function( 
 E_imports0
 , E_imports1
 , mu = 2021.299
 , sigma = 1/12/4
 , tfin = 2021.299+ 1.5/12
 , Rvariant = 2.0 
 , Rancestral = 1.5
 , rho0 = 0.1 
 , rho1 = 0.5 
) {
	library( lubridate )
	rho <- rho0 * rho1 
	nonext <- 5 / rho 
	
	itimes <- sim_importation_time(E_imports0, E_imports1, mu, sigma )
	
	## days
	.tfin <- as.numeric( as.Date( date_decimal( tfin ) ) - as.Date( '2021-01-01' ) )
	
	ctraj0 <- lapply( itimes$variant, function(tt){
		.tstart = as.numeric( as.Date( date_decimal( tt ) ) - as.Date( '2021-01-01' ) )
		sim_cluster_growth(.tstart, .tfin, Rvariant, non_extinction = 10 )
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
		sim_cluster_growth(.tstart, .tfin, Rancestral, non_extinction = 10 )
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


#~ o = sim_replicate( 
#~  E_imports0 = 25
#~  , E_imports1 = 25
#~  , mu = 2021.25
#~  , sigma = 1/12/4
#~  , tfin = 2021.25+ 2/12
#~  , Rvariant = 2.0 
#~  , Rancestral = 1.5
#~  , rho0 = 0.1 
#~  , rho1 = 0.5 
#~ )

#' Inference for simulated samples
#~ m = glm( (lineage=='variant') ~ sample_time , family = binomial( link = 'logit' ) , data = o ); summary( m )
