
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
sim_importation_time <- function(E_imports0, E_imports1, mu, sigma
	, shift_ancestral = 0#-21/365 
) 
{
	n0 <- rpois( 1, E_imports0 )
	n1 <- rpois( 1, E_imports1 )
	list(
		variant = rnorm( n0, mu, sigma )
		, ancestral = rnorm( n1, mu + shift_ancestral, sigma )
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
		cf = tryCatch( summary( m )$coefficients[ 2, 2 ], error = function(e){
			print(e)
			1
		})
		ms[[i]] <- c( 
			coef(m)[2] 
			, cf
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


#' @export 
sim_growthinference_hierbayes <- function( md, K = 10, iterations = 1e6, thin = 1e3 )
{
	library( glue )
	library( BayesianTools )
	vcluster <- md$cluster[ grepl(md$cluster, patt = '^variant' ) ]
	acluster <- md$cluster[ grepl(md$cluster, patt = '^ancestral' ) ]
	tvcluster <- table(vcluster)
	tacluster = table(acluster)
	keep <- c(
		names( tail( sort(tvcluster[ tvcluster > 1] ), K ))
		, names( tail( sort(tacluster[ tacluster > 1 ] ), K ))
		)
	md1 <- md[ md$cluster %in% keep , ]
	md1$icluster <- as.integer( md1$fcluster ) -1 
	
	## change timescale to generations 
	Tg <- 6.5
	md1$gtime <-  md1$sample_time / Tg 
	## make initial guess of alpha0, alpha1 and alpha2 
	Alpha01_0 <- t(sapply( keep, function( cid ){
		X = md1 
		X$y <- X$cluster == cid 
		m = glm( y ~ gtime, data = X , family = binomial( link = 'logit'))
		coef( m ) 
	}))
	
	## make initial parameter vector 
	#theta <- as.vector( Alpha01_0 )
	varnames <- keep[ grepl(keep, patt = '^variant') ]
	ancnames <- keep[ grepl(keep, patt = '^ancestral') ]
	alpha2_0 = unname( mean( Alpha01_0[varnames,2]  ) - mean( Alpha01_0[ancnames,2] ) ) 
	### recenter slope variables 
	Alpha01_0[ varnames, 2] <- Alpha01_0[ varnames, 2] - alpha2_0 
	
	keep1 <- keep[-1] # only estimate parms for k-1 
	
	intercept_names = paste(keep1, 'intercept', sep = '.')
	slope_names =  paste(keep1, 'slope', sep = '.')
	varnames1 <- intersect( varnames, keep1 )
	variant_slopes =  paste(varnames1, 'slope', sep = '.')
	theta <- c(  setNames( Alpha01_0[keep1,1], intercept_names  )
		, setNames( Alpha01_0[keep1,2], slope_names )
		, lineageEffect = alpha2_0 
	)
	M = length( theta )
	pnames <- names( theta )
	theta_upper = setNames( rep(max(2, theta[slope_names]), M), pnames )
	theta_upper[intercept_names] <- max( Inf, theta[intercept_names] )
	theta_lower = setNames( rep(min(-2, theta[slope_names]), M), pnames )
	theta_lower[intercept_names] <- min(-Inf, theta[intercept_names])
	taxis <- sort( unique( md1$gtime ))
	rprior <- function(n = 1)
	{
		if ( n > 1 ) {
			res =  matrix( setNames( rnorm( M, 0, 1 ), pnames ), ncol = 3, nrow = length(theta )  )
			colnames(res) <- names(theta)
		}
		setNames( rnorm( M, 0, 1 ), pnames )
	}
	dprior <- function(theta) {
		unname( dnorm( theta[ 'lineageEffect'] , 0, 1 , log = TRUE ) + 
		sum(dnorm( theta[slope_names] , 0, .1 , log = TRUE)) 
		)
	}
	dl <- function( .theta ) {
		
		theta1 <- theta
		theta1[ names(.theta) ] <- .theta 
		
		slope = theta1[ slope_names ]
		slope[variant_slopes] <- slope[variant_slopes] + theta1[ 'lineageEffect' ]
		inter = theta1[ intercept_names ]
		lox = sapply( 1:length( slope ), function(i){
			slope[i] * taxis + inter[i] 
		})
		# overdetermined
		px <- exp( lox ) /  ( 1 + exp( lox ))
		## fix any over 
		rspx <- rowSums(px)
		i <- which( rspx  > 1 ) 
		for ( ii in i ) 
			px[ii, ] <- px[ii, ] / rspx[ii] 
		## add in last cat 
		px <- pmax( px, 0 )
		px <- pmin( px, 1 )
		px <- cbind( px , 1 - rowSums( px ) )
		colnames( px  ) <- c( keep1 , keep[1] )
		coords = cbind(match(md1$gtime, taxis) , match( md1$cluster, colnames(px)) )
		pxterms = pmax(  px[ coords ] , 1e-4) # penalize zeros, but no inf's
		rv = sum( log ( pxterms )) #+ dprior( theta1 )
		#print(c( theta1, rv ))
		#print( unname( rv ) )
		if ( is.nan(rv))
			return( -Inf )
		rv
	}
	dl1 <- function( theta ){
		if ( is.matrix( theta )){
			return( apply( theta, MAR=2, dl ) )
		} 
		dl(theta )
	}
	# fit the model 
	pr0 <- createPrior( density = dprior, sampler = rprior, best = theta, upper = theta_upper, lower = theta_lower )
	bs0 <- createBayesianSetup(dl1, prior = pr0, names = names(theta)  ) #, parallel=3)
	.sv = matrix( rep( theta, 3 ) , ncol = 3 )
	rownames( .sv ) <- pnames 
	.sv <- t( .sv )
	f0 <- runMCMC(bs0, settings = list(iterations = iterations, thin = thin,  startValue = jitter(.sv))) 
	
	f0 
}

#' @export 
sim_growthinference_map <- function(md, K = 50 , bsrep = 100, ncpu = 1)
{
	library( glue )
	vcluster <- md$cluster[ grepl(md$cluster, patt = '^variant' ) ]
	acluster <- md$cluster[ grepl(md$cluster, patt = '^ancestral' ) ]
	tvcluster <- table(vcluster)
	tacluster = table(acluster)
	keep <- c(
		names( tail( sort(tvcluster[ tvcluster > 1] ), K ))
		, names( tail( sort(tacluster[ tacluster > 1 ] ), K ))
		)
	md1 <- md[ md$cluster %in% keep , ]
	
	## change timescale to generations 
	Tg <- 6.5
	md1$gtime <-  md1$sample_time / Tg 
	## make initial guess of alpha0, alpha1 and alpha2 
	Alpha01_0 <- t(sapply( keep, function( cid ){
		X = md1 
		X$y <- X$cluster == cid 
		m = glm( y ~ gtime, data = X , family = binomial( link = 'logit'))
		coef( m ) 
	}))
	
	## make initial parameter vector 
	varnames <- keep[ grepl(keep, patt = '^variant') ]
	ancnames <- keep[ grepl(keep, patt = '^ancestral') ]
	alpha2_0 = unname( mean( Alpha01_0[varnames,2]  ) - mean( Alpha01_0[ancnames,2] ) ) 
	### recenter slope variables 
	Alpha01_0[ varnames, 2] <- Alpha01_0[ varnames, 2] - alpha2_0 
	#Alpha01_0[ , 2] <- Alpha01_0[ , 2] - mean( Alpha01_0[ancnames,2] )  
	keep1 <- keep[-1] # only estimate parms for k-1 
	intercept_names = paste(keep1, 'intercept', sep = '.')
	slope_names =  paste(keep1, 'slope', sep = '.')
	varnames1 <- intersect( varnames, keep1 )
	variant_slopes =  paste(varnames1, 'slope', sep = '.')
	theta <- c(  setNames( Alpha01_0[keep1,1], intercept_names  )
		, setNames( Alpha01_0[keep1,2], slope_names )
		, lineageEffect = alpha2_0 
	)
	taxis <- sort( unique( md1$gtime ))
	dprior <- function(theta) {
		dnorm( theta[ 'lineageEffect'] , 0, 1 , log = TRUE ) + 
		sum(dnorm( theta[slope_names] , 0, .1 , log = TRUE)) 
	}
	of <- function( theta ) {
		slope = theta[ slope_names ]
		slope[variant_slopes] <- slope[variant_slopes] + theta[ 'lineageEffect' ]
		inter = theta[ intercept_names ]
		lox = sapply( 1:length( slope ), function(i){
			slope[i] * taxis + inter[i] 
		})
		# overdetermined
		px <- exp( lox ) /  ( 1 + exp( lox ))
		## fix any over 
		rspx <- rowSums(px)
		i <- which( rspx  > 1 ) 
		for ( ii in i ) 
			px[ii, ] <- px[ii, ] / rspx[ii] 
		## add in last cat 
		px <- pmax( px, 0 )
		px <- pmin( px, 1 )
		px <- cbind( px , 1 - rowSums( px ) )
		colnames( px  ) <- c( keep1 , keep[1] )
		coords = cbind(match(md1$gtime, taxis) , match( md1$cluster, colnames(px)) )
		pxterms = pmax(  px[ coords ] , 1e-4) # penalize zeros, but no inf's
		rv = sum( log ( pxterms )) + dprior( theta )
		#print(theta)
		#print( unname( rv ) )
		if ( is.nan(rv))
			return( -Inf )
		rv
	}
	o = optim( fn = of, par = theta, method = 'Nelder-Mead' , control = list(fnscale = -1, trace = 0, maxit = 5e3) , hessian=FALSE)
	
	bf <- c() 
	if ( bsrep > 0 ){
#~ browser()
		bf = .sim_inference_bsmap( md, nrep = bsrep, K = K, ncpu = ncpu  )
	}
	
	list( le = o$par['lineageEffect']
		, p = ifelse( length(bf)>0, mean( bf <= 0 ), 1 )
		, se = sd(bf) ,  keep = keep , fit = o$par , bs = bf )
}


.sim_inference_bsmap <- function ( md , K = 50, nrep = 500 , ncpu = 1)
{
	unlist(
	parallel::mclapply( 1:nrep , function(i){
		md1 <- md[ sample.int(nrow(md), replace=TRUE), ]
		f = sim_growthinference_map( md1, K = K, bsrep = 0, ncpu =1  )
		f$le 
	}, mc.cores = ncpu )
	)
}



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


#' Uses hier bayes model for identifying growth 
#'
#' @export 
sim_replicate2 <- function( MU = lubridate::decimal_date( as.Date( "2021-04-18")), TFIN, RHO0, k = 25, bs_replicates = 100,  ncpu = 1 , ... )
{
	o = sim_replicate( 
		  mu = MU #cf coguk/g2-adf3.rds, based on B117
		  , tfin = TFIN 
		  , rho0 = RHO0 
		 , ...
	)
	
	
	cdf <- o$clusterdf[ order( o$clusterdf$sample_time ) , ]
	cdf2 <- o$clusterdf2[ order( o$clusterdf$sample_time ) , ]
	
	
	f = tryCatch( sim_inference_clusterwise_logistic(o$clusterdf, minClusterSize = 5 , showres = FALSE)
	 , error = function(e) NULL )
	if ( is.null( f)){
		return(
			list( 
				data = cdf[1:min(nrow(cdf), 250e3) , ]
				, data2 = cdf[1:min(nrow(cdf2), 250e3) , ]
				, fit = data.frame( 
					tfin = TFIN 
					, window = floor( (TFIN - MU)*365 )
					, coef =0
					, p = 1
				)
				, fit_domestic = data.frame( 
					tfin = TFIN 
					, window = floor( (TFIN - MU)*365 )
					, coef = 0
					, p =1
				)
				, fit_bayes = data.frame( 
					tfin = TFIN 
					, window = floor( (TFIN - MU)*365 )
					, coef =0
					, p = 1
				)
				, mu = MU 
				, rho0 = RHO0 
				, minClusterSize = 5 
				, f1 = NA
				, f2 = NA 
				, X3 = NA
			)
		)
	}
	#c( (tt - MU)*365, o1$coef  )
	
	f2 = sim_inference_clusterwise_logistic(o$clusterdf2, minClusterSize = 5 , showres = FALSE)
	
	f3 = tryCatch( sim_growthinference_map(o$clusterdf, K = k , bsrep = bs_replicates, ncpu = ncpu )
	, error = function(e) {print(e); NULL} )
	if (is.null( f3 )){
		fbcoef <-  0
		fbp = 1
	} else{
		fbcoef <- f3$le
		fbp = f3$p
	}
	
	list( 
		data = cdf[1:min(nrow(cdf), 250e3) , ]
		, data2 = cdf[1:min(nrow(cdf2), 250e3) , ]
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
		, fit_bayes = data.frame( 
			tfin = TFIN 
			, window = floor( (TFIN - MU)*365 )
			, coef = fbcoef 
			, p = fbp
		)
		, mu = MU 
		, rho0 = RHO0 
		, minClusterSize = 5 
		, f1 = f 
		, f2 = f2 
		, f3 = f3 
		#, call = as.list(match.call())
	)
}


#' Three methods: gam adjusted for imports; simple logistic; hier bayes 
#' also simulates breakthrough 
#'
#' @export 
sim_replicate3 <- function( MU = lubridate::decimal_date( as.Date( "2021-04-18")), TFIN, RHO0, k = 25, bs_replicates = 100,  breakthroughrate = .41, or_breakthrough=1.25, ncpu = 1 , ... )
{
	o = sim_replicate( 
		  mu = MU #cf coguk/g2-adf3.rds, based on B117
		  , tfin = TFIN 
		  , rho0 = RHO0 
		 , ...
	)
	
	
	cdf <- o$clusterdf[ order( o$clusterdf$sample_time ) , ]
	cdf2 <- o$clusterdf2[ order( o$clusterdf$sample_time ) , ]
	
	
	## simulate breakthrough 
	or_breakthroughrate = breakthroughrate / ( 1- breakthroughrate )
	or_vbreakthroughrate = or_breakthrough * or_breakthroughrate 
	vbreakthroughrate <- or_vbreakthroughrate / ( 1 + or_vbreakthroughrate )
	
	cdf$breakthrough <- rbinom ( nrow(cdf) , size = 1, prob = breakthroughrate )
	i <- which( cdf$lineage == 'variant' )
	cdf[i,]$breakthrough <- rbinom ( length(i) , size = 1, prob = vbreakthroughrate )
	
	cdf2$breakthrough <- rbinom ( nrow(cdf2) , size = 1, prob = breakthroughrate )
	i <- which( cdf2$lineage == 'variant' )
	cdf2[i,]$breakthrough <- rbinom ( length(i) , size = 1, prob = vbreakthroughrate )
	
	### estimate or br
	brm = glm( breakthrough ~ lineage, data = cdf2 , family = binomial( link = 'logit' ))
	sbrm = summary( brm )
	
	
	## fit growth 
	library( mgcv )
	clusters <- unique( cdf$cluster )
	observed <- setNames( rep(FALSE, length( clusters )), clusters)
	cdf$first_obs <- FALSE
	cdf$first_obs_anc <- FALSE
	for ( i in 1:nrow( cdf ) ){
		if (   (cdf$lineage[i] == 'variant')  &    (!(observed[ cdf$cluster[i] ]))){
			cdf$first_obs[i] <- TRUE 
			observed[ cdf$cluster[i] ] <- TRUE 
		} else if ( (cdf$lineage[i] == 'ancestral')  &    (!(observed[ cdf$cluster[i] ])) ){
			cdf$first_obs_anc[i] <- TRUE 
			observed[ cdf$cluster[i] ] <- TRUE 
		}
	} 
	### variant importation intensity 
	f0 <- gam( first_obs ~ s(sample_time) , data = cdf[ cdf$lineage=='variant', ], family = binomial( link = 'logit'))
	f0a <- gam( first_obs_anc ~ s(sample_time) , data = cdf[ cdf$lineage=='ancestral', ], family = binomial( link = 'logit'))
	cdf$variant_importation_intensity = predict( f0, newdata = cdf  )
	cdf$ancestral_importation_intensity = predict( f0a, newdata = cdf  )
	cdf$net_import = with( cdf, exp( variant_importation_intensity  - ancestral_importation_intensity ) )
	f1 = glm ( (lineage=='variant') ~  net_import  + sample_time, family = binomial( link='logit'), data = cdf  )
	sf1 <- summary( f1 ) 
	
	f2 = glm( (lineage == 'variant') ~ sample_time , data = cdf , family = binomial( link = 'logit' ) )
	sf2 <- summary( f2 ) 
	
	f3 = tryCatch( sim_growthinference_map(o$clusterdf, K = k , bsrep = bs_replicates, ncpu = ncpu )
	, error = function(e) {print(e); NULL} )
	if (is.null( f3 )){
		fbcoef <-  0
		fbp = 1
	} else{
		fbcoef <- f3$le
		fbp = f3$p
	}
	
	list( 
		data = cdf[1:min(nrow(cdf), 250e3) , ]
		, data2 = cdf[1:min(nrow(cdf2), 250e3) , ]
		, fit_importation_adjusted = data.frame( 
			tfin = TFIN 
			, window = floor( (TFIN - MU)*365 )
			#~ 			, coef = unname( sf1$p.coef[ 'sample_time'])
			#~ 			, p = unname( sf1$p.pv[ 'sample_time' ] )
			, coef = unname( sf1$coefficients[3,1] ) 
			, p = unname( sf1$coefficients[3,4]  ) 
		)
		, fit_logistic = data.frame( 
			tfin = TFIN 
			, window = floor( (TFIN - MU)*365 )
			, coef = unname( sf2$coefficients[2,1] )
			, p = unname( sf2$coefficients[2,4]  ) 
		)
		, fit_bayes = data.frame( 
			tfin = TFIN 
			, window = floor( (TFIN - MU)*365 )
			, coef = fbcoef 
			, p = fbp
		)
		, fit_breakthrough = data.frame(
			tfin = TFIN 
			, window = floor( (TFIN - MU)*365 )
			, coef = unname( sbrm$coefficients[2,1] )
			, p = unname( sbrm$coefficients[2,4]  ) 
		)
		, mu = MU 
		, rho0 = RHO0 
		, minClusterSize = 5 
		#~ 		, f1 = f 
		#~ 		, f2 = f2 
		#~ 		, f3 = f3 
		#, call = as.list(match.call())
	)
}


