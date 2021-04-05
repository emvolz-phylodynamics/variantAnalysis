



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
sample_lineage <- function( root_dir = '..' ## '/cephfs/covid/bham/climb-covid19-volze'
	 , regenerate_ML_tree = F
	 , lineage = 'B.1.1.7' 
	 , weightsfn = '/cephfs/covid/bham/climb-covid19-volze/b0-weightsdf-2021-01-04.csv'
	 , mindate = as.Date( '2020-10-15' ) 
	 , maxdate = as.Date( Sys.Date() - 12 )
	 , nreps = 100 
	 , ns = c( 250, 500, 750, 1000 )
	 , prop_stratified = .25
	 , deduplicate = TRUE 
){
	library( lubridate )
	library( ape ) 
	library( glue )
	
	mintime <- decimal_date( mindate ) 
	maxtime <- decimal_date( maxdate )

	wdf = na.omit( read.csv( weightsfn , stringsAs=FALSE )  ) 


	civetfn =  list.files(  paste0(root_dir, '/phylolatest/civet/' ), patt = 'cog_global_[0-9\\-]+_metadata.csv', full.names=TRUE) #'../phylolatest/civet/cog_global_2020-12-01_metadata.csv'
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
	tr0 = read.tree( list.files( paste0(root_dir, '/phylolatest/trees'), patt = '.*newick', full=T ))
	tr0$tip.label <-  sapply( strsplit( tr0$tip.label, split='/' ), '[', 2 )
	tr1 <- keep.tip( tr0, s$central_sample_id )

	
	# t1 = Sys.time()
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
	# t2 = Sys.time()
	# t2-t1
	
	

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
	# n.b. Z[[1]] has length nreps
	
	if(regenerate_ML_tree) {
	  # algn_fn <- 	list.files(  paste0(root_dir, '/phylolatest/civet/' ), patt = 'cog_global_[0-9\\-]+_alignment.fasta', full.names=TRUE) #'../phylolatest/civet/cog_global_2020-12-01_metadata.csv'
	  # t0_readalgn = Sys.time()
	  # algn <- read.dna(algn_fn, format = 'fasta' ) 
	  # t1_readalgn = Sys.time()
	  # 
	  # t1_readalgn-t0_readalgn
	}
	
	
	
	for (X in Z ){
	  if(regenerate_ML_tree) {
  } else {
	    
	    Xs <- split( X, X$replicate )
	    tres <- lapply( Xs, function(y){
	      keep.tip( tr1, intersect(  y$central_sample_id, tr1$tip.label )  )
	    })
	    class( tres ) <- 'multiPhylo' 
	    write.tree( tres, file = glue('sampler1_{lineage}_{Sys.Date()}_n={X$sample_size[1]}{ifelse(deduplicate,"-deduped","")}.nwk')  )}
	}
	
	list(
	  trees = tres 
	  , tables = Z
	)
}



#' Matched sample by time (week) and UTLA
#'
#' @param csids central_sample_id to be matched 
#' @param lineage If not NULL, sample will only include this lineage 
#' @param not_lineage The sample will *not* include this lineage
#' @param nreps replicates
#' @param n_multiplier Optionally can draw multiple (integer) matches per element in csids 
#' @param deduplicate not implemented 
#' @export 
matched_sample <- function( 
  root_dir = '..' ## '/cephfs/covid/bham/climb-covid19-volze'
	, csids 
	 , lineage = NULL
	 , not_lineage = 'B.1.1.7'
	 , nreps = 1 
	 , n_multiplier = 1
	 , deduplicate = FALSE 
) {
  stopifnot( !deduplicate ) #not implemented

#~ csids <- read.tree('sampler1_B.1.1.7_2021-01-04_n=500-deduped.nwk')[[1]]$tip.label 
#~ lineage = NULL
#~ not_lineage = 'B.1.1.7'
#~ weightsfn = '/cephfs/covid/bham/climb-covid19-volze/b0-weightsdf-2021-01-04.csv'
#~ nreps = 100 
#~ n_multiplier = 1

	library( lubridate )
	library( ape ) 
	library( glue )
		
	civetfn =  list.files( paste0(root_dir, '/phylolatest/civet/') , patt = 'cog_global_[0-9\\-]+_metadata.csv', full.names=TRUE) #'../phylolatest/civet/cog_global_2020-12-01_metadata.csv'
	civmd = read.csv( civetfn , stringsAs=FALSE , header=TRUE )
	civmd$central_sample_id <-  sapply( strsplit( civmd$sequence_name , split='/' ) , '[', 2 ) # for linkage 
	civmd$sample_date <- as.Date( civmd$sample_date )
	civmd$sample_time <- decimal_date( civmd$sample_date ) 
	
	# load majora  '../latest/majora.20201204.metadata.matched.tsv' 
	jdf <- read.csv( list.files( paste0(root_dir, '/latest/' ), patt = 'majora.[0-9]+.metadata.matched.tsv', full.names=TRUE)  
	, stringsAs=FALSE, sep = '\t' ) 
	
	# combine
	s <- merge( civmd, jdf[ , c('central_sample_id', 'adm2') ], by = 'central_sample_id' )
	
	# exclude p1 
	lhls <- c( 'MILK', 'ALDP', 'QEUH', 'CAMC')
	s$lighthouse = FALSE
	for (lh in lhls ){
		s$lighthouse[ grepl(s$central_sample_id, patt = lh) ] <- TRUE
	}
	s <- s [ s$lighthouse , ]
	s$epiweek <- epiweek(  s$sample_date )
	s$year <- year ( s$sample_date )

	# load tree and deduplicate 
	tr0 = read.tree( list.files( paste0(root_dir,'/phylolatest/trees'), patt = '.*newick', full=T ))
	tr0$tip.label <-  sapply( strsplit( tr0$tip.label, split='/' ), '[', 2 )
	tr1 <- keep.tip( tr0, s$central_sample_id )
	
	#~ 	StatMatch package, NND.hotdeck
	## find controls 
	
	### filter lineage 
	s <- s[ , c('year' , 'epiweek', 'adm2', 'central_sample_id', 'lineage') ]
	s <- na.omit( s) 
	s$key <- with( s, paste( sep = '.' , year, epiweek, adm2 ) )
	srec <- s[ s$central_sample_id %in% csids , ] 
	sdon <- s[ !( s$central_sample_id %in% csids )  &  (s$lineage!=not_lineage) , ] 
	
	.sample <- function()
	{
		csid_control = do.call(c,  lapply( srec$key , function (key ){
			.s <- sdon[ sdon$key == key , ]
			n <- min( nrow( .s ) , n_multiplier ) 
			if ( n == 0 )
				return( NA )
			sample ( .s$central_sample_id, size = n , replace=FALSE )
		}) )
		csid_control <- na.omit( csid_control )
		csid_control
	}
	
	
	# sample over reps, ns; 
	{
		X = do.call( rbind, lapply(1:nreps, function(k){
			csids_control = .sample( )
			y = data.frame( central_sample_id = csids_control, replicate = k , sample_size = length( csids_control) )
			y
		}))
		write.csv( X, file = glue('{lineage}_matchSample_not{not_lineage}_{Sys.Date()}.csv')  )
	}
	#make tres
	{
		Xs <- split( X, X$replicate )
		tres <- lapply( Xs, function(y){
			keep.tip( tr1, intersect(  y$central_sample_id, tr1$tip.label )  )
		})
		class( tres ) <- 'multiPhylo' 
		write.tree( tres, file = glue('{lineage}_matchSample_not{not_lineage}_{Sys.Date()}.nwk')  )
	}
	
	list(
	 trees = tres 
	 , tables = X
	)
}

#' Matched sample by time (week) and UTLA for each sampled tree output from sample_lineage function
#'
#' @param sample_lineage_output output from sample_lineage
#' @param lineage If not NULL, will interpret as regular expression; only matches will be sampled
#' @param lineage_name incorporated into output file name 
#' @param not_lineage The sample will *not* include this lineage
#' @param n_multiplier Optionally can draw multiple (integer) matches per element in csids 
#' @param deduplicate not implemented 
#' @export 
matched_sample2 <- function( 
	sample_lineage_output_table
	 , lineage = NULL # 'B\\.1\\.177.*'  # pattern 
	 , lineage_name = NA  # for output file name
	 , not_lineage = 'B.1.1.7' # dont sample these
	 , n_multiplier = 1
	 , deduplicate = FALSE 
) {
stopifnot( !deduplicate ) #not implemented

	library( lubridate )
	library( ape ) 
	library( glue )
	
	csids_list =  lapply( split(sample_lineage_output_table, sample_lineage_output_table$replicate) , function(y) y$central_sample_id ) 
	nreps = length( csids_list )

	civetfn =  list.files(  '/cephfs/covid/bham/climb-covid19-volze/phylolatest/civet/' , patt = 'cog_global_[0-9\\-]+_metadata.csv', full.names=TRUE) #'../phylolatest/civet/cog_global_2020-12-01_metadata.csv'
	civmd = read.csv( civetfn , stringsAs=FALSE , header=TRUE )
	civmd$central_sample_id <-  sapply( strsplit( civmd$sequence_name , split='/' ) , '[', 2 ) # for linkage 
	civmd$sample_date <- as.Date( civmd$sample_date )
	civmd$sample_time <- decimal_date( civmd$sample_date ) 
	
	# load majora  '../latest/majora.20201204.metadata.matched.tsv' 
	jdf <- read.csv( list.files(  '/cephfs/covid/bham/climb-covid19-volze/latest/' , patt = 'majora.metadata.matched.tsv', full.names=TRUE)  
	, stringsAs=FALSE, sep = '\t' ) 
	
	# combine
	s <- merge( civmd, jdf[ , c('central_sample_id', 'adm2') ], by = 'central_sample_id' )
	
	# exclude p1 
	lhls <- c( 'MILK', 'ALDP', 'QEUH', 'CAMC')
	s$lighthouse = FALSE
	for (lh in lhls ){
		s$lighthouse[ grepl(s$central_sample_id, patt = lh) ] <- TRUE
	}
	s <- s [ s$lighthouse , ]
	s$epiweek <- epiweek(  s$sample_date )
	s$year <- year ( s$sample_date )

	# load tree and deduplicate 
	tr0 = read.tree( list.files( '/cephfs/covid/bham/climb-covid19-volze/phylolatest/trees', patt = '.*newick', full=T ))
	tr0$tip.label <-  sapply( strsplit( tr0$tip.label, split='/' ), '[', 2 )
	tr1 <- keep.tip( tr0, s$central_sample_id )
	
	#~ 	StatMatch package, NND.hotdeck
	## find controls 
	
	### filter lineage 
	s <- s[ , c('year' , 'epiweek', 'adm2', 'central_sample_id', 'lineage') ]
	s <- na.omit( s) 
	s$key <- with( s, paste( sep = '.' , year, epiweek, adm2 ) )
	
	
	.sample <- function(k)
	{
		csids = csids_list[[k]] 
		srec <- s[ s$central_sample_id %in% csids , ] 
		sdon <- s[ !( s$central_sample_id %in% csids )  &  (s$lineage!=not_lineage) , ] 
		if ( !is.null( lineage ) ){
			sdon <- sdon [ grepl(sdon$lineage, patt = lineage) , ]
		}
		
		csid_control = do.call(c,  lapply( srec$key , function (key ){
			.s <- sdon[ sdon$key == key , ]
			n <- min( nrow( .s ) , n_multiplier ) 
			if ( n == 0 )
				return( NA )
			sample ( .s$central_sample_id, size = n , replace=FALSE )
		}) )
		csid_control <- na.omit( csid_control )
		print(paste('made sample', k)) 
		csid_control
	}
	
	# sample over reps, ns; 
	{
		Xs <- lapply(1:nreps, function(k){
			csids_control = .sample( k )
			y = data.frame( central_sample_id = csids_control, replicate = k , sample_size = length( csids_control) )
			y
		})
		saveRDS( Xs, file = glue('matchSample_not{not_lineage}_lineage{lineage_name}_{Sys.Date()}.rds')  )
	}
	#make tres
	{
		tres <- lapply( Xs, function(y){
			keep.tip( tr1, intersect(  y$central_sample_id, tr1$tip.label )  )
		})
		class( tres ) <- 'multiPhylo' 
		write.tree( tres, file = glue('matchSample_not{not_lineage}_lineage{lineage_name}_{Sys.Date()}.nwk')  )
	}
	
	list(
	 trees = tres 
	 , tables_list = Xs
	)
}






