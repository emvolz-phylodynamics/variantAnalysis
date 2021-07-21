
#' Simulate cluster trajectories 
#'
#' SDE model for simple exponential growth, stores cumulative and active infections
#'
#' @param tstart Time of origin
#' @param tfin Simulate until 
#' @param R Expected reproduction number
#' @param Rsd std dev of R 
#' @param non_extinction Minimum cumulative infections required to return trajectory 
#' @param ... additional arguments are passed to seir model generator (e.g. for changing initial conditions)
#' @return Matrix with simulation output 
sim_cluster_growth <- function(tstart, tfin, R
	, Rsd = .2
	, non_extinction = 0
	, ntries = 10  
	, seir_gen = NULL 
	, ... )
{
	if ( tstart >= tfin ){
		return ( NULL )
	}
	.R <- max( 0, rnorm( 1, R , Rsd ))
	finci = -Inf 
	if ( is.null( seir_gen )){
		fn = system.file( 'stochastic_seir2.R', package = 'variantAnalysis' )
		suppressMessages( {seir_gen = odin::odin( fn ) } )
	}
	seir_sim <- seir_gen( R = .R ,  Requil = R, tini = tstart, tequil = tfin, ... )
	tries <- 0 
	while( (finci < non_extinction) & (tries < ntries) ) {
		X = seir_sim$run( seq( tstart, tfin, by = 1) ) 
		finci <- tail( X[, 'cI'], 1 )
		#print( finci )
		tries <- tries + 1
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
sim_sampling_cluster <- function( X,  rho0 = .10 , rho1 = 1-.66 ) 
{
	N <- floor( tail(X[, 'cI'], 1) )
	n <- rbinom( 1, N, rho0*rho1 )
	w <-  pmax(0, as.vector( X[, 'Il'] ) ) + pmax(0, as.vector( X[, 'Ih'] ) )
	lb <- min( 11, nrow(X)-1)
	taili <- (length(w ) - lb):length(w)
	w[taili] <- w[taili] *  exp( -(0:lb) / 3 )  # seq( 1-1/12, 0, by = -1/12 )
	if ( sum( w ) > 0 ){
		w  <- w / sum( w ) 
	} else{
		w <- rep(1 , length( w) )
	}
	if ( any( is.na( w ))) {
		#browser()  
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
#' @export 
sim_replicate <- function( 
 E_imports0 = 650 # cf coguk/g2-adf3.rds
 , E_imports1 = 650
 , mu = lubridate::decimal_date( as.Date( "2021-04-18")) #cf coguk/g2-adf3.rds, based on B117
 , sigma = 12.9 / 365 # std dev of cluster origin (yrs), g2-adf3.rds
 , tfin = lubridate::decimal_date( as.Date( "2021-04-18")) + 1/12
 , Rvariant = 1.5*1. #cf coguk/b.1.617/g3.1
 , Rancestral = 0.90 # approx for b117 in april
 , Rsd = .2 # std dev of initial R in clusters, cf coguk/b.1.617/g3.1
 , rho0 = 0.9 # proportion sequenced  
 , rho1 = 0.66 # proportion diagnosed 
 , seir_gen = NULL 
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
		sim_cluster_growth(.tstart, .tfin, Rvariant, Rsd= Rsd,  non_extinction = 10, seir_gen = seir_gen )
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
		ssc =  sim_sampling_cluster( ctraj1[[i]], rho0 = rho0, rho1 = rho1 ) 
		if ( length( ssc ) == 0 )
			return( NULL )
		data.frame( cluster = paste(sep='.', 'ancestral', i )
		 , sample_time = ssc
		 , lineage = 'ancestral'
		 , stringsAsFactors = FALSE
		)
	}))
	
	
	# also simulate large endogenous epidemic of incumbent lineage
	ninf0 <- length( itimes$variant ) * 5
	.tstart <- as.numeric( as.Date( date_decimal( min( itimes$variant)  ) ) - as.Date( '2021-01-01' ) ) - 14
	ctraj2 <- sim_cluster_growth(.tstart, .tfin, Rancestral, Rsd= 0.001, non_extinction = 10, EI0=ninf0 )
	ssc2 =  sim_sampling_cluster( ctraj2, rho0 = rho0, rho1 = rho1 ) 
	samp2 <- data.frame( cluster = paste(sep='.', 'ancestral', 'domestic' )
	 , sample_time = ssc2
	 , lineage = 'ancestral'
	 , stringsAsFactors = FALSE
	)
	list(
		clusterdf = rbind( samp0, samp1 ) # for comparison with reimported ancestral 
		, clusterdf2 = rbind( samp0, samp2 ) # for comparision with domestic ancestral 
		, variant_df = samp0
		, reimportedAncestral_df = samp1
		, domesticAncestral_df = samp2 
	)
}






#~ #' Inference for simulated samples

#~ library( mlogit ) 
#~ minsize <- 20 
#~ tc <- table( o$cluster )
#~ keepc <- names(tc) [ tc >= minsize ] 
#~ o1 <- o[ o$cluster %in% keepc , ]

#~ m = glm( (lineage=='variant') ~ sample_time , family = binomial( link = 'logit' ) , data = o ); summary( m )
#~ m = mgcv::gam( as.factor( cluster )  ~ sample_time * lineage , data = o1 , family = mgcv::multinom(K = length(unique( o1$cluster)))  )


#' @export 
sim_inference_clusterwise_logistic <- function(s, minClusterSize = 25 , showres = FALSE)
{
	# filter by cluster size 
	x = table( s$cluster ) ; x = x[ x >= minClusterSize ]; keep<- names(x) 
	s <- s[ s$cluster %in% keep , ]
	
	# compute spans and make sure they cover period ; exclude others 
	s_clusts <- split( s, s$cluster )
	clusts <- names( s_clusts )
	clustspans <- sapply( s_clusts , function(ss){
		range( ss$sample_time )
	})# time range of cluster 
	colnames( clustspans ) <- clusts 
	
	clust_lineage <- sapply( s_clusts, function(x)  x$lineage[1]  )
	
	
	s_genotype = split( s, s$lineage )
	
	# analysis 1, individual cluster level 
	ms <- list() # fits 
	for(i in 1:length(clusts))
	{
	  cc <- clusts[i]
		ss = s_clusts[[ cc ]] 
		ccg = ss$lineage[1] 
		notccg = setdiff( c('ancestral', 'variant') , ccg )
		X <- rbind( ss , s_genotype[[ notccg ]] )
		X <- X[ with(X, sample_time>=clustspans[1,cc] & sample_time <=clustspans[2,cc]) ,]
		X$y <- X$lineage == ccg 
		
		m = glm( y ~ sample_time , data = X, family = binomial(link='logit' ))
		
		# store to plot
		ms[[i]] <- c( 
			coef(m)[2] 
			, summary( m ) $coef[ 2,2] 
		)
	}
	rs_ests <- do.call( cbind , ms ) 
	rs_ests_wt <- rs_ests[ , clust_lineage == 'ancestral' ]
	rs_ests_mutant <- rs_ests[ , clust_lineage == 'variant' ]
	
	if ( showres ){
		library( sinaplot )
		sinaplot( list( rs_ests_wt[1,], rs_ests_mutant[1, ] ) )
	}
	#print( kruskal.test( list(  rs_ests_wt[1,], rs_ests_mutant[1, ] )  )  )
	
	d <- as.data.frame( t( rs_ests ))
	colnames(d) <- c( 's', 'stderr'  )
	d$lineage <- clust_lineage
	d$w <- with( d , 1 / stderr )
	
	m1 = lm( s ~ lineage , data =d, weight = w )
	summ1 = summary( m1 )
	list(
		coef = summ1$coefficients[ 2, c(2, 4)]
		, fitcoefs = rs_ests 
		, res = d
		#, summary = summ1 
	)
}
#~ sim_inference_clusterwise_logistic( o  )



#' @export 
sim_replicate1 <- function( MU = lubridate::decimal_date( as.Date( "2021-04-18")), TFIN, RHO0, ... )
{
	o = sim_replicate( 
		  mu = MU #cf coguk/g2-adf3.rds, based on B117
		  , tfin = TFIN 
		  , rho0 = RHO0 
		 , ...
	)
	f = sim_inference_clusterwise_logistic(o$clusterdf, minClusterSize = 5 , showres = FALSE)
	#c( (tt - MU)*365, o1$coef  )
	
	f2 = sim_inference_clusterwise_logistic(o$clusterdf2, minClusterSize = 5 , showres = FALSE)
	
	print( nrow(o) )
	list( 
		data = o$clusterdf[1:min(nrow(o$clusterdf), 100e3) , ]
		, fit = data.frame( 
			tfin = TFIN 
			, window = floor( (TFIN - MU)*365 )
			, coef = unname( f$coef[1] )
			, p = unname( f$coef[2]  ) 
		)
		, fit_domestic = data.frame( 
			tfin = TFIN 
			, window = floor( (TFIN - MU)*365 )
			, coef = unname( f2$coef[1] )
			, p = unname( f2$coef[2]  ) 
		)
		, mu = MU 
		, rho0 = RHO0 
		, minClusterSize = 5 
		, f1 = f 
		, f2 = f2 
		#, call = as.list(match.call())
	)
}


