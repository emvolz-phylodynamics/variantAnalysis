library( BayesianTools ) 

#' Estimate rate that clusters increase or decrease in frequency and transform into selection coefficients for different genotypes
#' 
#' @param s data frame with sequence_name(character), del_introduction(character), sample_time(numeric), genotype (wt or mutant)
#' @param mint lower sample time bound. Clusters must have sample before this time
#' @param maxt upper sample time bound. Clusters must have sample after this time
#' @param minClusterSize exclude clusters with fewer samples than this
#' @param r specifies generation time (1/r) when transforming to selection coefficients. Default 73 corresponds to five days.
#' @return  list with matrix of selection coefficients for different clusters and a t.test for difference in 'mutant' selection coefficients
#' @export 
clusterwise_logistic <- function(s, mint = 2020.20, maxt = 2020.35, minClusterSize = 25 , r = 73 , quasi=F)
{
	# filter missing
	s <- s[with( s, !is.na(sample_time) & !is.na(del_introduction) & !is.na(genotype) ) , ]
	
	# filter by cluster size 
	x = table( s$del_introduction ) ; x = x[ x >= minClusterSize ]; keep<- names(x) 
	s <- s[ s$del_introduction %in% keep , ]
	
	# compute spans and make sure they cover period ; exclude others 
	s_clusts <- split( s, s$del_introduction )
	clusts <- names( s_clusts )
	clustspans <- sapply( s_clusts , function(ss){
		range( ss$sample_time )
	})# time range of cluster 
	colnames( clustspans ) <- clusts 
	keep <- colnames( clustspans )[ clustspans[1,] <= mint  &  clustspans[2,] >= maxt ]
	s <- s[ s$del_introduction %in% keep ,  ]
	
	# remove sample times outside of mint and maxt 
	s <- s[ s$sample_time >= mint   , ] #&  s$sample_time <= maxt
	if ( nrow(s) < 100 ) 
		stop('mint and maxt too restrictive' )
	s <- s[ , colnames(s)%in%c('del_introduction', 'sample_time', 'genotype','weight') ] ## remove identifiers 
	rownames(s) <- NULL 
	# recompute clust spans 
	s_clusts <- split( s, s$del_introduction )
	clusts <- names( s_clusts )
	clustspans <- sapply( s_clusts , function(ss){
		range( ss$sample_time )
	})# time range of cluster 
	colnames( clustspans ) <- clusts 
	
	# majority rule for genotype
	clust_genotypes <- sapply( s_clusts, function(x) {
		tt = table(x$genotype)
		names( tt[ which.max(tt) ] )[1] 
	})
	
	s_genotype = split( s, s$genotype )
	
	# analysis 1, individual cluster level 
	ms <- list()
	for(i in 1:length(clusts)){
	  cc <- clusts[i]
		ss = s_clusts[[ cc ]] 
		ccg = ss$genotype[1] 
		notccg = setdiff( c('wt', 'mutant') , ccg )
		X <- rbind( ss , s_genotype[[ notccg ]] )
		X <- X[ with(X, sample_time>=clustspans[1,cc] & sample_time <=clustspans[2,cc]) ,]
		X$y <- X$genotype == ccg 
		if(quasi==F){
		  m = glm( y ~ sample_time , data = X, family = binomial(link='logit' ))
		}else{
		  m = glm( y ~ sample_time , data = X, family = quasibinomial(link='logit' ),weight=weight)
		}
		# store to plot
		ms[[i]] <- m
	}
	rs_ests <- sapply(ms,function(m) c( coef(m)[2] , tryCatch({confint(m)[2,]}, error = function(e) { return(c(NA,NA))})  ))
	rs_ests_wt <- rs_ests[ , clust_genotypes == 'wt' ]
	rs_ests_mutant <- rs_ests[ , clust_genotypes == 'mutant' ]
	
	#sinaplot::sinaplot( list( rs_ests_wt[1,], rs_ests_mutant[1, ] ) )
	#kruskal.test( list(  rs_ests_wt[1,], rs_ests_mutant[1, ] )  ) 
	s_mutant = rs_ests_mutant/r 
	smut = s_mutant[1, ] 
	w = apply( s_mutant , 2, function(x) 1/abs( diff( tail(x,2)) )) 
	m = lm( smut ~ 1, weight = w ); 
	
	list( s_mutant = s_mutant
		, s_wt = rs_ests_wt/r
		, selcoef = c( coef(m) , confint( m )  )
		, t.test = t.test ( rs_ests_mutant[1, ]/r  ) 
		, data = s
		, data_clusts  = s_clusts 
		, clust_genotypes = clust_genotypes 
		,	ms = ms
	)
}



#' Hierarchical Bayesian inference of cluster-wise selection coefficients with a random effect corresponding to a dichotomous genotype
#' 
#' @param s data frame that describes each sample, must have columns: del_introduction, sample_time, genotype 
#' @param MINT minimum sample time for including cluster 
#' @param MAXT exclude samples after this date 
#' @param MU sets time scale of process and generation time. Generation time is 1/MU years
#' @param minClusterSize minimum size for cluster to be included in analysis 
#' @return Bayesian sampler output produced by BayesianTools::runMCMC
#' @export 
hier_bayes_exponentialGrowth_frequency <- function(s,  iterations = 8000000, thin = 200 , ncpu = 4, MINT =  -Inf, MAXT = 2020.25, MU = 73, minClusterSize = 10)
{
	
	library( BayesianTools )
	s <- s[ s$sample_time > MINT  &  s$sample_time <= MAXT  , ] 
	s <- s[ order( s$sample_time ) , ]
	x = table( s$del_introduction ) ; x = x[ x >= minClusterSize]; keep_lineages <- names(x) 
	s <- s [ s$del_introduction %in% keep_lineages , ]
	
	#  make d:  data frame that describes each cluster, must have columns: lineage (matching del_introduction), first_sample (time of first sample in cluster), genotype ( 'wt' or 'mutant' )
	Xs <- split( s, s$del_introduction )
	# majority rule for genotype
	clust_genotypes <- sapply( Xs, function(x) {
		tt = table(x$genotype)
		names( tt[ which.max(tt) ] )[1] 
	})
	d = data.frame(
		lineage = sapply( Xs, function(x) x$del_introduction[1] )
		,genotype =  clust_genotypes
		,first_sample  = sapply( Xs, function(x) min( x$sample_time) )
		, stringsAsFactors=FALSE
	)
	
	# remove any clusters without samples 
	d <- d[ d$lineage %in% unique( s$del_introduction ) , ] 
	M <- nrow(d) 
	# define index of each lineage 
	d$ilineage <- 1:nrow(d)
	s$ilineage <- d$ilineage[ match( s$del_introduction, d$lineage )]

	# indicator for which lineages are extant at each sample 
	EXTANT <- matrix (0, nrow = nrow(s), ncol = nrow(d) )
	for ( k in 1:nrow(d)){
		EXTANT[ s$sample_time >= d$first_sample[k] , k ] <- 1
	}
	# t minus start of cluster 
	TminusTi <- matrix (-Inf, nrow = nrow(s), ncol = nrow(d) )
	for ( k in 1:nrow(d)){
		TminusTi[ ,k] <-  s$sample_time  - d$first_sample[k] 
	}
	
	# prior 
	## initial parameters and bounds 
	pnames <- c( paste(sep='.', 'logn', 1:M) , paste( sep='.', 'r', 1:M), 'cv' )
	theta0 <- c( rep(0, M)
	  , rep( 0, M )
	  ,  1
	)
	theta0.ub <- c( rep(2, M)
	  , rep( MU, M )
	  ,  5
	)
	theta0.lb <- c( rep(-2, M)
	  , rep(-MU, M )
	  ,  .01
	)
	names(theta0) = names(theta0.lb) = names(theta0.ub) = pnames
	## density
	dprior0 <- function(theta){
		theta <-  pmax( theta0.lb, pmin( theta0.ub, theta ))
		logns <- theta[1:M] 
		rs <- theta[ (M+1):(M+M) ]
		theta1 <- tail( theta, 1 ) 
		cv <- theta1[1]
		sigma <- cv * MU 
		
		drs <- rs[ d$genotype=='wt' ]
		grs <- rs[ d$genotype=='mutant' ]
		
		dmu = mean( drs ) 
		gmu = mean( grs ) 
		
		.selcoef <- ( gmu - dmu ) / MU 
		
		pr0 <- sum( dnorm( drs , dmu, sigma, log=TRUE ) ) + sum( dnorm(grs, gmu, sigma, log=TRUE) )
		
		pr1 <- sum( dnorm( logns, 0, 1 , log=TRUE ) )
		
		pr2 <- dnorm( .selcoef, 0, 1, log =TRUE) +  dgamma(cv , shape=1, scale=1 , log=TRUE)
		
		pr3 <- dnorm( gmu / MU, 0, 1, log=TRUE) + dnorm( dmu/MU, 0, 1, log=TRUE )
		
		unname( pr0 + pr1 + pr2 + pr3)
	}
	## sampler for starting conditions 
	rprior0 <- function(n=1) 
	{
		if ( n > 1){
			return( 
				sapply( 1:n, function(k) 
				  pmax( theta0.lb, pmin( theta0.ub,  
						 c( rnorm(M,0,.25)
						  , rnorm( M, 0, MU/5)
						  , runif(1, .5, 2)
						)
				  )) 
				 ) 
			)
		}
	  pmax( theta0.lb, pmin( theta0.ub,  
			 c( rnorm(M,0,.25)
			  , rnorm( M, 0, MU/5 )
			  , runif(1, .5, 2)
			)
	  )) 
	}

	# objective functions 
	ll0 <- function( theta ) 
	{
		theta <-  pmax( theta0.lb, pmin( theta0.ub, theta ))
		logns <- theta[1:M] 
		rs <- theta[ (M+1):(M+M) ]
		
		denmat <- t(EXTANT) * exp( logns ) * exp(t(TminusTi) * rs )
		den <- log( colSums( denmat ) )
		
		num <- ( logns[ s$ilineage ]  +   TminusTi[  cbind(1:nrow(s), s$ilineage)   ] * rs[s$ilineage] )
		ll <- sum(num) - sum(den) 
		ll
	}
	ll0.1 <- function( theta ){
		if ( is.matrix( theta )){
			return( apply( theta, MAR=2, ll0 ) )
		} 
		ll0(theta )
	}

	
	# fit the model 
	pr0 <- createPrior( density = dprior0, sampler = rprior0,  lower = theta0.lb, upper = theta0.ub, best = theta0 )
	bs0 <- createBayesianSetup(ll0.1, prior = pr0, names = pnames, parallel=ncpu)
	f0 <- runMCMC(bs0, settings = list(iterations = iterations, thin = thin,  startValues = bs0$prior$sampler()))
	i <- Sys.getpid() 
	saveRDS( f0, file=paste0('a5f0_', i, '.rds' ) )
	
	f0$clusterdata <- d 
	f0$data = s
	f0$clust_genotypes = clust_genotypes
	f0 
}






#' too slow 
hier_bayes_pairwiseLogistic_frequency <- function(s, d
 , iterations = 8000000, thin = 200 
, ncpu = 4
, r = 73 # 5 days 
, minclustsize = 10 
)
{

	MINT <- -Inf 
	MAXT <- 2020.25

	s <- s[ s$sample_time > MINT  &  s$sample_time <= MAXT  , ] 
	s <- s[ order( s$sample_time ) , ]
	x = table( s$del_introduction ) ; x = x[ x >= minclustsize]; keep_lineages <- names(x) 
	s <- s [ s$del_introduction %in% keep_lineages , ]

	# remove any clusters without samples 
	d <- d[ d$lineage %in% unique( s$del_introduction ) , ] 
	M <- nrow(d) 
	# define index of each lineage 
	d$ilineage <- 1:nrow(d)
	s$ilineage <- d$ilineage[ match( s$del_introduction, d$lineage )]
	
	nc = nrow(d) 
	n = nrow(s) 

	# indicator for which lineages are extant at each sample 
	EXTANT <- matrix (0, nrow = nrow(s), ncol = nrow(d) )
	for ( k in 1:nrow(d)){
		EXTANT[ s$sample_time >= d$first_sample[k] , k ] <- 1
	}
	# t minus start of cluster ( lineages X cluster )
	TminusTi <- matrix (-Inf, nrow = nrow(s), ncol = nrow(d) )
	for ( k in 1:nrow(d)){
		TminusTi[ ,k] <-  s$sample_time  - d$first_sample[k] 
	}
	ti.c.c <- array( NA, c( nrow(d), nrow(d), nrow(s)  ))
	for ( i in 1:n ) 
	{
		x = TminusTi[i, ]
		ti.c.c[,,i] <- sapply( 1:nc, function(k) sapply( 1:nc, function(l) min( x[k], x[l] )  ))
	}
	
	# compute sel coef between clusters 
	.smat <- function ( fs ) {
		fs %*% t(1/fs) - 1
	}
	
	# prior 
	dprior0 <- function(theta){
		theta <-  pmax( theta0.lb, pmin( theta0.ub, theta ))
		p0s <- theta[1:M] 
		fs <- theta[ (M+1):(M+M) ]
		theta1 <- tail( theta, 1 ) 
		logsd <- theta1[1]
		
		wtfs = fs[ d$genotype == 'wt' ]
		mutfs = fs [ d$genotype == 'mutant'] 
		wtmu = mean( wtfs )
		mutmu = mean( mutfs )
		
		.selcoef <- mutmu / wtmu - 1
		
		pr0 <- sum( dnorm( log(wtfs), log(wtmu), logsd, log=TRUE ) ) +  sum( dnorm( log(mutfs), log(mutmu), logsd, log=TRUE ) )
		
		pr1 <- sum( dbeta( p0s, 1, 4 , log=TRUE ) )
		
		pr2 <- dnorm( .selcoef, 0, 1, log =TRUE) +  dexp( logsd , rate = 1/5,  log=TRUE )
		
		unname( pr0 + pr1 + pr2 )
	}
	
	# initial parameters and bounds 
	pnames <- c( paste(sep='.', 'p0', 1:M) , paste( sep='.', 'f', 1:M), 'logsd' )
	theta0 <- c( rep(1/M, M)
	  , rep(1, M )
	  ,  .5
	)
	theta0.lb <- c( rep(1e-6, M)
	  , rep(.1, M )
	  ,  .01
	)
	theta0.ub <- c( rep(1, M)
	  , rep( 10, M )
	  , 1
	)
	names(theta0) = names(theta0.lb) = names(theta0.ub) = pnames

	# sampler for starting conditions 
	rprior0 <- function(n=1) 
	{
		if ( n > 1){
			return( 
				sapply( 1:n, function(k) 
				  pmax( theta0.lb, pmin( theta0.ub,  
						 c( rbeta(M,1,4)
						  , rlnorm( M, 0, .25 )
						  , runif(1, .25, .75)
						)
				  )) 
				 ) 
			)
		}
		pmax( theta0.lb, pmin( theta0.ub,  
			 c( rbeta(M,1,4)
						  , rlnorm( M, 0, .25 )
						  , runif(1, .25, .75)
						)
		)) 
	}
	
	
	# objective functions 
	ll0 <- function( theta ) 
	{
		theta <-  pmax( theta0.lb, pmin( theta0.ub, theta ))
		p0s <- theta[1:M] 
		fs <- theta[ (M+1):(M+M) ]
		smat <- .smat ( fs )
		
		x0kl <- matrix( 1, nrow = nc , ncol = nc ) 
		nextant <- 1
		isextant <- rep(FALSE, nc )
		ki = s$ilineage[1]
		isextant[ki] <- TRUE 
		p <- rep(0, nc ); p[ki] <- p0s[ki]
		ll = 0
		for( i in 2:nrow(s)){
			ki = s$ilineage[i]
			ti = s$sample_time[i] 
			iextant <- which(isextant)
			if ( nextant > 1 ){
				A <- matrix( 0, ncol = nc , nrow = choose( nextant, 2 ) +1 )
				A[1, iextant] <- 1 
				#apply(  t( combn(iextant,2) ) , 1,  function(kl) {
				#	k = kl[1]
				#	l = kl[2] 
				ai = 2
				for ( .k in 1:(length( iextant )-1)) for ( .l in (.k+1):length(iextant)) {
					k <- iextant[.k]
					l <- iextant[.l]
					Arow = rep(0, ncol (A) ) 
					Arow[k] <- 1
					Arow[l] <- -x0kl[k,l] * exp( r * smat[k,l] * ti.c.c[ k, l, i]  )
					A[ai,] <- Arow 
					ai <- ai + 1 
				}
				b <- rep(0, nrow(A)); b[1] <- 1
				p = coef( lm( b ~ A-1 )  ) 
				p [is.na(p)] <- 0 
				p <- p / sum(p) 
			}
			
			p0i = p[ki] 
			if ( ti == d$first_sample[ki] ) {
				p0i <- p0s[ ki ]
				p <- (1-p0i) * p /  sum(p) 
				p[ki] <- p0i
				nextant <- nextant +1 
				isextant [ ki] <- TRUE 
				for ( .l in 1:(length( iextant )))  {
					l <- iextant[.l]
					x0kl[ki, l] <- p0i / p[l] 
					x0kl[ l, ki] <- p[l] / p0i 
				}
#~ browser() 
			}
			if ( nextant > 1 ){
				ll <- ll + log( p0i) - log( 1-p0i)
			}
		}
		ll
	}

#~ ll0( rprior0() )
#~ stop() 

	ll0.1 <- function( theta ){
		if ( is.matrix( theta )){
			return( apply( theta, MAR=2, ll0 ) )
		} 
		ll0(theta )
	}



	# fit the model 
	pr0 <- createPrior( density = dprior0, sampler = rprior0,  lower = theta0.lb, upper = theta0.ub, best = theta0 )
	bs0 <- createBayesianSetup(ll0.1, prior = pr0, names = pnames, parallel=ncpu)
	f0 <- runMCMC(bs0, settings = list(iterations = iterations, thin = thin,  startValues = bs0$prior$sampler()))
	i <- Sys.getpid() 
	saveRDS( f0, file=paste0('hier_bayes_pairwiseLogistic_frequency', '_', i, '.rds' ) )
# TODO include d in output, indicator for genotype 
	f0 
}


#' Compute cluster sizes only counting clusters with first sample occuring within given range and arrange by genotype 
#' 
#' @param s data frame with sequence_name(character), del_introduction(character), sample_time(numeric), genotype (wt or mutant)
#' @param mint lower sample time bound. Clusters must have sample before this time
#' @param maxt upper sample time bound. Clusters must have sample after this time
#' @return data frame with cluster sizes, origin time, and genotype 
#' @export 
cluster_sizes <- function(s, mint = 2020.20, maxt = 2020.35)
{
#~  mint = decimal_date( as.Date('2020-08-01')) 
#~  maxt = decimal_date( as.Date('2020-10-07'))
	# filter missing
	s <- s[with( s, !is.na(sample_time) & !is.na(del_introduction) & !is.na(genotype) ) , ]
	
	# compute spans and make sure they **begin** within period ; exclude others 
	s_clusts <- split( s, s$del_introduction )
	clusts <- names( s_clusts )
	clust_starts <- sapply( s_clusts , function(ss){
		min( ss$sample_time )
	})
	names( clust_starts ) <- clusts 
	keep <- names( clust_starts )[ clust_starts >= mint  &  clust_starts <= maxt ]
	s <- s[ s$del_introduction %in% keep ,  ]
	
	# remove sample times outside of mint and maxt 
	s <- s[ s$sample_time >= mint  &  s$sample_time <= maxt , ]
	s <- s[ , c('del_introduction', 'sample_time', 'genotype') ] ## remove identifiers 
	rownames(s) <- NULL 
	# recompute clust spans 
	s_clusts <- split( s, s$del_introduction )
	clusts <- names( s_clusts )
	clustspans <- sapply( s_clusts , function(ss){
		range( ss$sample_time )
	})# time range of cluster 
	colnames( clustspans ) <- clusts 
	clust_starts <- sapply( s_clusts , function(ss){
		min( ss$sample_time )
	})
	names( clust_starts ) <- clusts 
	clust_genotypes <- sapply( s_clusts, function(x) { 
		# majority rule 
		tx = table( x$genotype )
		names(tx)[which.max(tx)]
	})
	
	s_genotype = split( s, s$genotype )
	
	X= data.frame( cluster = clusts 
	 , origin = clust_starts
	 , size = sapply( s_clusts, nrow )
	 , genotype = clust_genotypes 
	)
	X
}

#' @export 
plot_cluster_sizes <- function(X, mincs = 2 ){
	library( ggplot2 )
	X1 <- X [ X$size  > mincs , ]
	p = ggplot( X1 , aes( x = origin, y = size, colour = genotype ) ) + geom_point()  + scale_y_log10() + stat_smooth( method = lm ) 
	p
}



#' Performs a simple logistic model fit and make a frequency plot 
#' 
#' genotype has two values: 'ancestral'(or 'wt') and mutant 
#'
#' @param s1 data frame with columns genotype and sample_time 
#' @param cols Colours for plotting 
#' @param bw bin width for aggregating data 
#' @export 
plpg <- function( s1, bw = 1 , cols = c(ancestral='white', mutant='black')
		 , g = 365 / 6.5 # recovery rate 
		 , RD = seq( 1.1,1.5, length = 100) #hypothetical rep numbers for ancestral 
		 , gentime = 6.5/365 
		 , timeScaleMethod = c( 'gentime', 'growth') 
		)
{ 
	library( ggplot2 )
	library( lubridate ) 
	library( cowplot )
	
	x = table( s1$genotype  )
	if ( 'wt' %in% names (x ) )
		s1$genotype [ s1$genotype == 'wt'] <- 'ancestral'
	
	m = glm( (genotype=='mutant') ~ sample_time, family = binomial(link='logit'),  data = s1  )
	print( summary( m ))
	
	timeScaleMethod = timeScaleMethod[1] 
	rs = c( coef( m )[2] ,  confint( m ) [2, ]  )
	if (timeScaleMethod=='growth') {
		rD <- RD * g - g
		#~ selcoef <- rs / rD 
		selcoef <- (rs) %*% t( 1/rD )
		print(median( selcoef ))
		print(range( selcoef ))
		selcoef = c( median( selcoef ), range( selcoef ) )
	} else if ( timeScaleMethod=='gentime' ){
		selcoef = rs * gentime 
	}
	
	s1 <- s1[ !is.na( s1$sample_time ), ]
	.cols = cols
	.time <- seq( min( s1$sample_time )-bw/365, max( s1$sample_time )+bw/365, by = bw / 365 )
	h = hist( s1$sample_time,  breaks = .time , plot=FALSE)
	s1dg = split( s1, s1$genotype )
	hd = hist( s1dg$ancestral$sample_time, breaks = .time , plot=F)
	hg = hist( s1dg$mutant$sample_time, breaks = .time , plot=F)
	
	sweek = data.frame( time = hd$mids, nD = hd$counts, nG = hg$counts, date = date_decimal( hd$mids )  )
	sweek$n = with( sweek, nD + nG )
	sweek$pG <- with( sweek, nG / n )
	
	
	#.time <- seq( 2020, 2020.25, length = 100 )
	theta2pg <- function( theta ){
		t50 = -theta[1] / theta[2]
		logodds50 = theta[1]+t50 * theta[2]
		#pg50 = exp( logodds50 ) / ( 1 + exp( logodds50) ) # = 50%, checking the math 
		pg.time <- exp( theta[2] * (sweek$time - t50 ) ) / ( 1 +  exp( theta[2] * (sweek$time - t50 ) ) )
		pg.time
	}
	sweek$fitted <- theta2pg( coef(m) )
	
	prof = profile( m, which =1  )[[1]]
	thetaub = tail( unlist( tail( prof[ prof$z < -1.96 , ], 1) ), 2 )
	thetalb = tail( unlist( head( prof[ prof$z  > 1.96 , ], 1) ), 2 )
	sweek$fittedub = theta2pg( thetaub )
	sweek$fittedlb = theta2pg( thetalb )
	
	sweek$nglb =  with( sweek , qbinom( .025, size = n, prob = fittedub)  )
	sweek$ngub =  with( sweek, qbinom( .975, size = n, prob = fittedlb )  )
	sweek$pglb <- with( sweek, nglb / n )
	sweek$pgub <- with( sweek, ngub / n )
	k = with( sweek ,  pmin(nD,nG)==0 )
	sweek$pgub[ k ] <- 1
	sweek$pglb[ k ] <- 0
	
	p0 = ggplot( data = sweek, aes( x = date, y = fitted , ymax = pgub, ymin=pglb ) ) + geom_path() + geom_point( aes( x = date, y = pG, size =n ), colour = 'black' ) + geom_ribbon(alpha = .2 ) + theme_classic()  + theme( legend.position = 'none' ) + xlab('') + ylab('Frequency of sampling mutant')
	
	X = with (sweek, data.frame( date = date, n = nD , genotype='ancestral' ) )
	X <- rbind( X, with (sweek, data.frame( date = date, n = nG , genotype='mutant' ) ) )
	p1 = ggplot( data= X , aes( x = date, y = n, fill = genotype)) + geom_bar( position = 'stack', stat='identity', colour = 'grey' , alpha=.5) + scale_fill_manual( values = .cols)  + xlab('') + ylab('Samples') + theme_classic() + theme( legend.position='top')
	p = plot_grid( plotlist= list( p0 , p1 ), nrow =  2) # ggtitle( paste('sel coef', selcoef))
	p$m = m 
	list( selcoef = selcoef , plot1=p0, plot2=p1, plot3=p, modelfit=m )
}



#' Compute cumulative proportion of samples descended from clusters originating within a specified time window
#'
#' @param s data frame with  del_introduction(character), sample_time(numeric)
#' @param mint lower sample time bound. Clusters must have sample before this time
#' @param maxt upper sample time bound. Clusters must have sample after this time
#' @param res integer resolution of time axis
#' @return list with 1) result with  columns: cumulative fraction of samples descended from new clusters, number descended from old clusters, number descended from new clusters 2) data with cluster subset 
#' @export
cluster_origin_comparison <- function(s, mint = decimal_date( as.Date('2020-08-01')) , maxt = Inf , res = 50 )
{
	if ( 'del_lineage' %in% colnames(s) & !('del_introduction' %in% colnames(s) ) )
			s$del_introduction <- s$del_lineage
	# filter missing
	s <- s[with( s, !is.na(sample_time) & !is.na(del_introduction) ) , ]

	# compute spans and make sure they cover period ; exclude others
	s_clusts <- split( s, s$del_introduction )
	clusts <- names( s_clusts )
	clustspans <- sapply( s_clusts , function(ss){
			range( ss$sample_time )
	})# time range of cluster
	colnames( clustspans ) <- clusts
	newClusts <-  clusts[ clustspans[1,] >= mint & clustspans[2, ] <= maxt ]
	s$isnew <- s$del_introduction %in% newClusts

	# remove sample times outside of mint and maxt
	s <- s[ s$sample_time >= mint  &  s$sample_time <= maxt , ]
	rownames(s) <- NULL
	s <- s[ order(s$sample_time ) , ]

	ss <- split( s, s$isnew )

	taxis <- seq( mint, max( s$sample_time ), length = res )
	X = t( sapply( taxis , function(tt){
			c( sum( ss[[1]]$sample_time <= tt ), sum(  ss[[2]]$sample_time <= tt ) )
	}) )
	colnames(X) <- names( ss )
	pnew = X[,'TRUE'] / rowSums( X )


	res = cbind( taxis, pnew, X )
	list( result = res, data = s )
}



#' Compute cumulative proportion of samples descended from clusters originating within a specified time window or uk lineage 
#' 
#' This returns a data frame with lieange log odds broken down by epi week 
#'
#' @param s data frame with  del_introduction(character), sample_time(numeric), and epi_week (integer)
#' @param uk_lineage  optional character vector of uk_lineage
#' @param mint lower sample time bound. Clusters must have sample before this time
#' @param maxt upper sample time bound. Clusters must have sample after this time
#' @return data frame with log odds of new lineages or uk_lineage 
#' @export
cluster_origin_comparison2 <- function(s, uk_lineage = NULL, genotype = NULL, mint = decimal_date( as.Date('2020-07-01')) , maxt = Inf, detailed=TRUE  )
{
	if ( 'del_lineage' %in% colnames(s) & !('del_introduction' %in% colnames(s) ) )
			s$del_introduction <- s$del_lineage
	# filter missing
	s <- s[with( s, !is.na(sample_time) & !is.na(del_introduction) ) , ]

	# compute spans and make sure they cover period ; exclude others
	s_clusts <- split( s, s$del_introduction )
	clusts <- names( s_clusts )
	clustspans <- sapply( s_clusts , function(ss){
			range( ss$sample_time )
	})# time range of cluster
	colnames( clustspans ) <- clusts
	newClusts <-  clusts[ clustspans[1,] >= mint & clustspans[2, ] <= maxt ]
	s$isnew <- s$del_introduction %in% newClusts

	# remove sample times outside of mint and maxt
	s <- s[ s$sample_time >= mint  &  s$sample_time <= maxt , ]
	rownames(s) <- NULL
	s <- s[ order(s$sample_time ) , ]
	
	if (!is.null( uk_lineage) )
		ss = split( s, s$uk_lineage==uk_lineage  )
	else if ( !is.null( genotype ))
		ss = split( s, s$genotype == genotype  )
	else
		ss <- split( s, s$isnew )
	
	weeks = seq( min ( s$epi_week ) , max( s$epi_week ))
#~ browser() 
	weights = sapply( weeks, function(w){
		n1 <- sum( ss[['TRUE']]$epi_week == w  )   
		n2 = sum( ss[['FALSE']]$epi_week == w  )  
		n <- n1 + n2 
		if ( n1 == 0 ) 
			return(0)
		if ( n2 == 0 )
			return(0) 
		p = sum( ss[['TRUE']]$epi_week == w  )  / n
		v = p * (1 - p) / n 
		1 / sqrt( v )
	})
	weights = weights / median( weights )
	logodds = sapply( weeks, function(w){
		n1 <- sum( ss[['TRUE']]$epi_week == w  )   
		n2 = sum( ss[['FALSE']]$epi_week == w  )  
		if ( n1 == 0 | n2 == 0 )
			return(NA)
		log( n1 )  - log(n2 )
	})
	
	
	#list( result=cbind(  weeks = weeks, weights = weights, logodds = logodds ) , data = s )
	res = data.frame(  weeks = weeks, weights = weights, logodds = logodds )
	if ( detailed) 
		res = list(
		  result = res 
		  , data = s
		)
	return(res)
}




#' Simple log odds freqency plot for given variable. Does not partition by cluster. 
#'
#' If the column 'weight' is in the input data frame, this will use the values to weight each observation. 
#' If computing for genotypes you will want to remove gaps and 'X' before passing to this function; also deduplicate by patient id. 
#' 
#' @param s data frame which must contain sample_time and 'variable' for each sequence
#' @param variable Column of s which contains variable to split the data 
#' @param value The value in conlumn 'variable' for which we will compute frequency 
#' @param mint Minimum sample time to include
#' @param maxt Max time to include
#' @param form formula passed to mgcv::gam for estimating frequencies 
#' @return list with data, plots, and estimated frequencies 
#' @export 
variable_frequency_epiweek <- function(s, variable='genotype', value='mutant',  mint = -Inf , maxt = Inf, detailed=TRUE , form = y~s(sample_time, bs = 'gp', k = 5) )
{
	library( lubridate ) 
	# remove sample times outside of mint and maxt
	if ( ( 'sample_date' %in% colnames(s) )  & !('sample_time' %in% colnames(s) ) )
	{
		s$sample_time <- decimal_date( as.Date( ymd( s$sample_date ) )  )
	}
	stopifnot( is.numeric( s$sample_time ) )
	s <- s[ s$sample_time >= mint  &  s$sample_time <= maxt , ]
	rownames(s) <- NULL
	s <- s[ order(s$sample_time ) , ]
	
	if ( !('epi_week' %in% colnames(s)))
		s$epi_week <- lubridate::epiweek( s$sample_time )
	
	weeks = seq( min ( s$epi_week ) , max( s$epi_week ))
	weekstarts = sapply( weeks, function(x) min( na.omit( s$sample_time[ s$epi_week==x ] ) ) )
	weekends = sapply( weeks, function(x) max( na.omit( s$sample_time[ s$epi_week==x ] ) ) )
	
	if ( 'weight' %in% colnames(s))
		s$weight <- s$weight / mean( s$weight ) 
	
	ss <- split( s, s[[variable]]==value )

	weights = sapply( weeks, function(w){
		if ( 'weight' %in% colnames(s) ) {
			n1 <- sum( ss[['TRUE']]$weight[ ss[['TRUE']]$epi_week==w ]  )
			n2 <- sum( ss[['FALSE']]$weight[ ss[['FALSE']]$epi_week==w ]  )
		} else{
			n1 <- sum( ss[['TRUE']]$epi_week == w  )   
			n2 = sum( ss[['FALSE']]$epi_week == w  )  
		}
		n <- n1 + n2 
		if ( n1 == 0 ) 
			return(0)
		if ( n2 == 0 )
			return(0) 
		p = sum( ss[['TRUE']]$epi_week == w  )  / n
		v = p * (1 - p) / n 
		1 / sqrt( v )
	})
	weights = weights / median( weights )
	logodds = sapply( weeks, function(w){
		if ( 'weight' %in% colnames(s) ) {
			n1 <- sum( ss[['TRUE']]$weight[ ss[['TRUE']]$epi_week==w ]  )
			n2 <- sum( ss[['FALSE']]$weight[ ss[['FALSE']]$epi_week==w ]  )
		} else{
			n1 <- sum( ss[['TRUE']]$epi_week == w  )   
			n2 = sum( ss[['FALSE']]$epi_week == w  )  
		}
		if ( n1 == 0 | n2 == 0 )
			return(NA)
		log( n1 )  - log(n2 )
	})
	
	library( mgcv ) 
	.s1 <- s[ order( s$sample_time ) , ]
	.s1 <- .s1[ .s1$sample_time >= (min(ss[['TRUE']]$sample_time)-1/365), ]
	.s1$y <- .s1[[variable]]==value
	m = mgcv::gam( form , family = binomial(link='logit') , data = .s1 )
	.s1$estimated = predict( m )
	estdf = data.frame( time = .s1$sample_time, estimated_frequency = .s1$estimated )
	
	res = data.frame(  time = weekends, weeks = weeks, weights = weights, logodds = logodds )
	if ( detailed) {
		library( ggplot2 )
		res = list(
		  result = res 
		  , estimated =  estdf
		  , data = .s1
		  , plot =  ggplot(  aes( x = as.Date( date_decimal( time ) ), y = logodds, size=weights, weight=weights ) , data= res ) + geom_point ()  + theme_minimal() + theme( legend.pos='') + geom_smooth(method=stats::loess) + xlab('') + ggtitle( paste0('Frequency of ', variable, '=', value) ) #, method.args=list(span = 1) 
		)
	}
	return(res)
}
