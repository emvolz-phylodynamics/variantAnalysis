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
clusterwise_logistic <- function(s, mint = 2020.20, maxt = 2020.35, minClusterSize = 25 , r = 73 )
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
	clust_genotypes <- sapply( s_clusts, function(x) x$genotype[1] )
	
	s_genotype = split( s, s$genotype )
	
	# analysis 1, individual cluster level 
	rs_ests <- sapply( clusts , function(cc){
		ss = s_clusts[[ cc ]] 
		ccg = ss$genotype[1] 
		notccg = setdiff( c('wt', 'mutant') , ccg )
		X <- rbind( ss , s_genotype[[ notccg ]] )
		X <- X[ with(X, sample_time>=clustspans[1,cc] & sample_time <=clustspans[2,cc]) ,]
		X$y <- X$genotype == ccg 
		m = glm( y ~ sample_time , data = X, family = binomial(link='logit' ))
		c( coef(m)[2] , confint( m ) [2, ]  )
	})
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
	)
}



#' Hierarchical Bayesian inference of cluster-wise selection coefficients with a random effect corresponding to a dichotomous genotype
#' 
#' @param s data frame that describes each sample, must have columns: del_introduction, sample_time
#' @param d data frame that describes each cluster, must have columns: lineage (matching del_introduction), first_sample (time of first sample in cluster), genotype ( 'wt' or 'mutant' )
#' @param MINT minimum sample time for including cluster 
#' @param MAXT exclude samples after this date 
#' @param MU sets time scale of process and generation time. Generation time is 1/MU years
#' @return Bayesian sampler output produced by BayesianTools::runMCMC
#' @export 
hier_bayes_exponentialGrowth_frequency <- function(s, d, iterations = 8000000, thin = 200 , ncpu = 4, MINT =  -Inf, MAXT = 2020.25, MU = 73)
{
	library( BayesianTools )
	s <- s[ s$sample_time > MINT  &  s$sample_time <= MAXT  , ] 
	s <- s[ order( s$sample_time ) , ]
	x = table( s$del_introduction ) ; x = x[ x >= 10]; keep_lineages <- names(x) 
	s <- s [ s$del_introduction %in% keep_lineages , ]

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
	  , rep( 10, M )
	  ,  5
	)
	theta0.lb <- c( rep(-2, M)
	  , rep(-10, M )
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
		
		.selcoef <- ( dmu - gmu ) / MU 
		
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
	f0 
}
