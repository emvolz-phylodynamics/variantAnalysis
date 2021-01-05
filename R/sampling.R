



#' Produce a random sample of given global lineage and sample weights
#'
#' Sampling is repeated multiple times and over multiple sample sizes. Can confine sampling to a particular time frame. A proportion of samples can be stratified through time rather than being drawn completely at random (weighted smapling is still used within each week)
#'
#' @param lineage Global lineage
#' @param weightsfn output produced by coverage_weights
#' @param mindate minimum date
#' @param maxdate maximum date; note it takes about 12 days for P2 sampling to stabilize 
#' @param nreps replicates
#' @param ns sample sizes to use 
#' @param prop_stratified proportion of sample to reserve for stratified sampling 
#' @param deduplicate if TRUE, duplicate samples will be removed 
#' @export 
sample_lineage <- function( 
	lineage = 'B.1.1.7' 
	 , weightsfn = '/cephfs/covid/bham/climb-covid19-volze/b0-weightsdf-2021-01-04.csv'
	 , mindate = as.Date( '2020-10-15' ) 
	 , maxdate = as.Date( Sys.Date() - 12 )
	 , nreps = 100 
	 , ns = c( 250, 500, 750, 1000 )
	 , prop_stratified = .25
	 , deduplicate = TRUE 
) {
	library( lubridate )
	library( ape ) 
	library( glue )
	
	mintime <- decimal_date( mindate ) 
	maxtime <- decimal_date( maxdate )

	wdf = na.omit( read.csv( weightsfn , stringsAs=FALSE )  ) 


	civetfn =  list.files(  '../phylolatest/civet/' , patt = 'cog_global_[0-9\\-]+_metadata.csv', full.names=TRUE) #'../phylolatest/civet/cog_global_2020-12-01_metadata.csv'
	civmd = read.csv( civetfn , stringsAs=FALSE , header=TRUE )
	civmd$central_sample_id <-  sapply( strsplit( civmd$sequence_name , split='/' ) , '[', 2 ) # for linkage 
	civmd$sample_date <- as.Date( civmd$sample_date )
	civmd$sample_time <- decimal_date( civmd$sample_date ) 

	# filter dates 
	civmd <- civmd[ civmd$sample_time >= mintime , ]
	civmd <- civmd[ civmd$sample_time <= maxtime , ]

	# filter lineage 
	civmd <- civmd[ civmd$lineage == lineage , ]

	# combine
	s <- merge( civmd, wdf, by = 'central_sample_id' )

	# exclude p1 
	lhls <- c( 'MILK', 'ALDP', 'QEUH', 'CAMC')
	s$lighthouse = FALSE
	for (lh in lhls ){
		s$lighthouse[ grepl(s$central_sample_id, patt = lh) ] <- TRUE
	}
	s <- s [ s$lighthouse , ]
	s$epiweek <- epiweek(  s$sample_date )

	# load tree and deduplicate 
	tr0 = read.tree( list.files( '../phylolatest/trees', patt = '.*newick', full=T ))
	tr0$tip.label <-  sapply( strsplit( tr0$tip.label, split='/' ), '[', 2 )
	tr1 <- keep.tip( tr0, s$central_sample_id )

	if ( deduplicate )
	{
		if ( Ntip (tr1) > 6e3 ){
			stop('Sample size is large. Current dedup code cannot handle this size: ', Ntip(tr1 )) 
		}
		D <- as.matrix( cophenetic.phylo( tr1 )  ) #  TODO not sustainable 
		diag(D) <- 0 
		clusters = sapply(rownames(D), function(x) { 
			paste( sort( colnames(D)[ D[x,] < 1e-6 ] ), collapse = '.' )
		})
		s$cluster = NA 
		s$cluster <- clusters[ match( rownames(D), s$central_sample_id ) ]
		s <- s[ order( s$sample_time ) , ]
		.s <- s[ !duplicated( s$cluster ), ]
	}

	# main sampling func
	.sample <- function( n ) 
	{
		nstrat <- floor ( n * prop_stratified )
		nn <- n - nstrat
		npw <- floor( nstrat / length( unique( s$epiweek )) ) 
		## make strat sample
		stratsample <- do.call( c, lapply( split( s, s$epiweek ), function(ss){
			sample( ss$central_sample_id , prob = ss$coverage_weight, replace=FALSE, size = min(nrow(ss),npw) )
		}))
		## make remainder 
		sids <- setdiff( s$central_sample_id , stratsample )
		s1 <- s[ !( s$central_sample_id %in% stratsample ), ]
		osample <- with( s1, sample( central_sample_id , prob = coverage_weight, replace=FALSE, size = min( nrow(s1), n - length(stratsample)) ))
		c( stratsample, osample )
	}
	#~ .sample( 100 )


	# sample over reps, ns; 
	Z <- list() 
	for ( n in ns ){
		X = do.call( rbind, lapply(1:nreps, function(k){
			y = data.frame( central_sample_id = .sample( n ), replicate = k , sample_size = n )
		}))
		Z[[ as.character(n) ]] <- X
		write.csv( X, file = glue('sampler1_{lineage}_{Sys.Date()}_n={X$sample_size[1]}{ifelse(deduplicate,"-deduped","")}.csv')  )
	}
	#make tres
	for (X in Z ){
		Xs <- split( X, X$replicate )
		tres <- lapply( Xs, function(y){
			keep.tip( tr1, intersect(  y$central_sample_id, tr1$tip.label )  )
		})
		class( tres ) <- 'multiPhylo' 
		write.tree( tres, file = glue('sampler1_{lineage}_{Sys.Date()}_n={X$sample_size[1]}{ifelse(deduplicate,"-deduped","")}.nwk')  )
	}
	
	list(
	 trees = tres 
	 , tables = Z
	)
}
