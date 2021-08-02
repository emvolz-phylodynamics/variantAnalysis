#' Compute scanner statistics for nodes in tree. 
#'
#' Takes standard inputs in the form of a rooted phylogeny and data frame with required metadata (see below). 
#' Computes a logistic growth rate and a statistic for outlying values in a molecular clock (root to tip regression).
#' Comparison sample is matched by time and region.
#' 
#' add optional 'sample_weight' to metadata, use in regression models 
#' 
#' @param tre A rooted phylogeny in ape::phylo or treeio::treedata form. The latter may include internal node labels and be read by '.read.beast'. 
#' @param amd A data frame containing required metadata for each tip in tree: sequence_name, sample_date, region. Optional metadata includes: pillar_2, sample_time(numeric), country. 
#' @param min_descendants Clade must have at least this many tips 
#' @param max_descendants Clade can have at most this many tips 
#' @param min_date Only include samples after this data 
#' @param max_date Only include samples before and including this date
#' @param min_cluster_age_yrs Only include clades that have sample tips that span at least this value
#' @param ncpu number cpu for multicore ops 
#' @param output_dir Path to directory where results will be saved 
#' @param include_pillar1 if TRUE (default TRUE), will include Pillar 1 samples when computing stats
#' @param country Optional, only include data from this country 
#' @param num_ancestor_comparison When finding comparison sample for molecular clock stat, make sure sister clade has at least this many tips 
#' @param factor_geo_comparison  When finding comparison sample based on geography, make sure sample has this factor times the number within clade of interest
#' @param Tg Approximate generation time in years 
#' @param report_freq Print progress for every n'th node 
#' @export 
scanner2 <- function(tre
 , amd
 , min_descendants = 30 
 , max_descendants = 20e3
 , min_cluster_age_yrs = 1.5/12
 , min_date = NULL, max_date = NULL , ncpu = 8
 , output_dir = '.' 
 , include_pillar1 = TRUE 
 , country = NULL 
 , num_ancestor_comparison = 500
 , factor_geo_comparison = 5
 , Tg = 6.5/365
 , report_freq = 50 
)
{
print(paste('Starting ', Sys.time()) )
	#library( treeio ) 
	library( ape ) 
	library( lubridate )
	library( glue ) 
	
	max_time <- Inf 
	if ( !is.null( max_date )){
		max_time <- decimal_date( max_date )
	} else{
		max_date <- Sys.Date()
	}
	
	min_time <- -Inf 
	if (!is.null( min_date ))
		min_time <- decimal_date( min_date )
	
	# tree data 
	if ( inherits ( tre, 'treedataList' ) | inherits( tre, 'treedata' )){
		nodedata = as.data.frame( get.data( tre ) )
		tre = get.tree( tre )
		rownames( nodedata ) <- nodedata$node
	} else if ( inherits( tre, 'phylo' )){
		nodedata = NULL 
	} else{
			stop('tre must be of type ape::phylo or treeio::treedata')
	}
	# load coguk algn md 
	amd <- amd[ !is.na( amd$sequence_name ) , ]
	amd$sample_date <- as.Date( amd$sample_date )
	if ( !('sample_time' %in% colnames(amd))){
		amd$sample_time = decimal_date (amd$sample_date)
	}
	amd$sts <- amd$sample_time 
	# exclude missing dates 
	amd <- amd[ !is.na( amd$sample_time ) , ] 
	if ( !('country' %in% colnames(amd)) )
		amd$country = 'country_not_specified'
	if ( !('pillar_2' %in% colnames(amd)) ){
	  amd$pillar_2 = 'True'
	  if( ('is_pillar_2' %in% colnames(amd)) ){ 
	    amd$pillar_2 = amd$is_pillar_2
	    amd$pillar_2 = ifelse(amd$is_pillar_2 == "Y", "True", ifelse(amd$is_pillar_2 == "N", "False", NA))
	  } 
	  
	}
	# only data within country 
	if ( !is.null( country ))
		amd <- amd[ amd$country == country, ]
	# exclude p1 
	if ( !include_pillar1 )
	{
		amd <- amd[ amd$pillar_2 == 'True' , ]
	}
	
	# filter by sample time 
	amd <- amd [ (amd$sample_time >= min_time) & (amd$sts <= max_time) , ] 
	sts <- setNames( amd$sample_time , amd$sequence_name )
	
	## treat tip labs not in amd differently 
	tre$tip.label [ !(tre$tip.label %in% amd$sequence_name) ] <- NA 
	
	# root to tip 
	ndel <- node.depth.edgelength( tre ) 
	
	message(paste('Loaded tree & filtered by inclusion criteria', Sys.time()) )
	
	# data structures to quickly look up tree data 
	# copied/adapted from treestructure 
	{
		n <- Ntip( tre ) 
		nnode = Nnode( tre )
		poedges <- tre$edge[ ape::postorder( tre ), ]
		preedges <- tre$edge[ rev( ape::postorder( tre )), ]
		
		
		# PRECOMPUTE for each node 
		# num descendents
		ndesc <- rep( 0, n + nnode ) # n tips descending 
		ndesc[1:n] <- 1 
		Ndesc <- rep( 1, n + nnode ) # include internal nodes
		# max and min sample times of descendants :
		tlsts <- sts[ tre$tip.label ]
		max_desc_time <- rep(-Inf, n + nnode ); max_desc_time[1:n] <- tlsts 
		min_desc_time <- rep(Inf, n + nnode ); min_desc_time[1:n] <- tlsts 
		# postorder traverse
		for (ie in 1:nrow(poedges)){
			a <- poedges[ie,1]
			u <- poedges[ie,2]
			ndesc[a] <- ndesc[a] + ndesc[u] 
			Ndesc[a] <- Ndesc[a] + Ndesc[u]  
			max_desc_time[a] <- max( max_desc_time[a] , max_desc_time[u] )
			min_desc_time[a] <- min( min_desc_time[a] , min_desc_time[u] )
		}
		clade_age <- max_desc_time - min_desc_time # time span of descendant tips
		
		
		# descendants; note also counts mrca node 
		descendants <- lapply( 1:(n+nnode), function(u) integer( Ndesc[u] ) )  #pre allocate 
		for ( u in 1:(n+nnode) )
			descendants[[u]][1] <- u 
		Ndesc_index <- rep( 2, n + nnode )
		for (ie in 1:nrow(poedges)){
			a <- poedges[ie,1]
			u <- poedges[ie,2]
			## mem efficient way to fill in values for a
			i0 <- Ndesc_index[a] 
			i1 <- Ndesc[u] + i0 - 1
			descendants[[a]][ i0:i1 ] <- descendants[[u]]  
			Ndesc_index[a] <- i1 + 1 
		}
		
		# descendant tips; index of tip 
		descendantTips <- descendants
		for (a in (n+1):(n+nnode))
			descendantTips[[a]] <- descendantTips[[a]][ descendantTips[[a]] <= n ]
		## tip labels, only including names in amd which meet inclusion criteria (excl NA )
		descendantSids <- lapply( 1:(n+nnode), function(u) na.omit(  tre$tip.label[ descendantTips[[u]] ] ) )
		
		# ancestors
		st0 <- Sys.time()
		ancestors <- lapply( 1:(n+nnode), function(u) integer() )
		for (ie in 1:nrow(preedges)){
			a <- preedges[ie,1]
			u <- preedges[ie,2]
			ancestors[[u]] <- c( ancestors[[a]], a)
		}
		st1 <- Sys.time()
		
		## vector samp time desc 
		descsts = NULL #deprecate 
		
	}
	message(paste('Derived lookup variables', Sys.time()) ) 
	
	
	.get_comparator_ancestor <- function(u, num_comparison = num_ancestor_comparison)
	{
		nu = ndesc[u] 
		asu = ancestors[[u]]
		na = -1
		for ( a in asu ){
			na = ndesc[a]
			nc = na - nu 
			if ( nc >= num_comparison )
			{
				break 
			}
		}
		if ( na < nu ){
			if ( na > 0 ){
				message('Failed to find comparator ancestor with more tips. Returning NA. Node: ', u) 
			}
			return (NA) 
		}
		a
	}
	
	# matched by time and in proportion to region prevalence
	.get_comparator_sample <- function( u , nX = factor_geo_comparison ) 
	{
		# weight for region 
		 w = table ( amd$region[ match( descendantSids[[u]]  , amd$sequence_name ) ]  )  
		 w = w [ names(w)!='' ]
		 w = w / sum( w ) 
		 
		 nu <- ndesc[u] 
		 tu =  descendantSids[[u]] 
		 stu = sts[ tu ]
		#~ 		 stu = descsts [[ u ]]
		 minstu = min(na.omit(stu ))
		 maxstu = max(na.omit(stu) )
		 
		 amd1 = amd[ (amd$sample_time >= minstu) & (amd$sample_time <= maxstu), c('sequence_name', 'sample_date', 'region') ]
		 amd1 <- amd1[ !(amd1$sequence_name %in% tu) , ]
		 amd1 <- amd1[ amd1$region %in% names( w ) , ]
		 na = min ( nrow( amd1 ) , nu * nX )
		 if ( na < nu ) 
			return ( NULL )
		 amd1$w = w[ amd1$region ] 
		 ta = sample( amd1$sequence_name, replace=FALSE, size = na , prob = amd1$w ) 
		 ta
	}
	
	.logistic_growth_stat <- function(u, ta = NULL, generation_time_scale = Tg)
	{
		if ( is.null(ta))
			ta = .get_comparator_sample(u) 
		if ( is.null( ta ))
			return(0)
		tu = descendantSids[[u]] 
		sta = sts[ ta ]
		stu = sts[ tu ]
		X = data.frame( time = c( sta, stu ), type = c( rep('control', length(ta)), rep('clade',length(tu)) ) )
		X = na.omit( X ) 
		m = glm( type=='clade' ~ time, data = X, family = binomial( link = 'logit' ))
		s = summary( m ) 
		rv = unname( coef( m )[2] * generation_time_scale ) 
		p = NA 
		if ( is.na( rv )){
			message( 'NA growth stat, node: ', u  )	
		} else{ 
			p = s$coefficients[2, 4 ]
		}
		c( rv, p )
	}
	
	# log median p-value of rtt predicted divergence of tips under u
	.clock_outlier_stat <- function(u)
	{
		a = .get_comparator_ancestor(u)
		if ( is.na( a ))
			return( NA ) 
		
		tu = descendantSids[[u]]
		ta = setdiff( descendantSids[[a]] , tu )
		
		sta = sts[ ta ]
		stu = sts[ tu ]
		#~ 		stu = descsts [[ u ]]
		
		iu = match( tu , tre$tip.label )
		ia = match( ta , tre$tip.label )
		ndelu = ndel[ iu ] 
		ndela = ndel[ ia ] 
		m = lm ( ndela ~ sta ) 
		r2 = summary( m )$r.squared
		oosp = predict(m, newdata =  data.frame(sta = unname( stu )) )  -  ndelu
		#mean( oosp) / sqrt( mean( (predict(m) - ndela )^2 ) )
		sqrt( mean( oosp^2 )  )
	}
	
	.lineage_summary <- function(tips, maxrows = 4){
		if ( is.null(tips))
			return( '' )
		lins = amd$lineage[ match( tips, amd$sequence_name)]
		tx = sort(table( lins), decreasing=TRUE ) / length( lins ) 
		if ( length( tx ) >1 ){
			y = as.data.frame( tx ) 
			colnames(y) <- c( 'Lineage', 'Frequency')
		} else{
			y = data.frame( Lineage = lins[1], Frequency = 1 )
		}
		y <- y [ 1:min(nrow(y),maxrows) , ]
		y$Frequency <- paste0( round(y$Frequency*100), '%' )
		paste( knitr::kable(y, 'simple') , collapse = '\n' ) # convert to string
	}
	
	.region_summary <- function(tips, maxrows = 5){
		if ( is.null(tips))
			return( '' )
		regs = amd$region[ match( tips, amd$sequence_name)]
		tx = sort(table( regs ), decreasing=TRUE ) / length( regs ) 
		if ( length( tx ) >1 ){
			y = as.data.frame( tx ) 
			colnames(y) <- c( 'Region', 'Frequency')
		} else{
			y = data.frame( Region = regs[1], Frequency = 1 )
		}
		y <- y [ 1:min(nrow(y),maxrows) , ]
		y$Frequency <- paste0( round(y$Frequency*100), '%' )
		paste( knitr::kable(y, 'simple') , collapse = '\n' ) # convert to string
	}
	
	# main 
	## compute stats for subset of nodes based on size and age 
	nodes = which(  (ndesc >= min_descendants)   &   (ndesc <= max_descendants) & (clade_age >= min_cluster_age_yrs) )
	report_nodes <- nodes[ seq(1, length(nodes), by = report_freq) ]
	Y = do.call( rbind, 
		parallel::mclapply( nodes , function(u){
			tu = descendantSids[[u]]
			ta = .get_comparator_sample(u) 
			ulins <- amd$lineage[ match( tu, amd$sequence_name)]
			alins <- amd$lineage[ match( ta, amd$sequence_name)]
			lgs = .logistic_growth_stat ( u, ta )
			X = data.frame( cluster_id = ifelse(is.null(nodedata)
				, as.character(u) 
				,  as.character(nodedata[ as.character(u), 'cluster_id' ]) 
				) 
			 , node_number = u 
			 , parent_number = ifelse( is.null(ancestors[[u]] ), NA, tail( ancestors[[u]], 1 ) ) 
			 , most_recent_tip = as.Date( date_decimal( max( na.omit(sts[ tu ])  ) ) )
			 , least_recent_tip = as.Date( date_decimal( min( na.omit( sts[ tu ])  ) ) )
			 , cluster_size = length( tu )
			 , logistic_growth_rate = lgs[1]
			 , logistic_growth_rate_p = lgs[2] 
			 , clock_outlier = .clock_outlier_stat(u)
			 , lineage = paste( names(sort(table(ulins),decreasing=TRUE)) , collapse = '|' )
			 , lineage_summary = .lineage_summary( tu ) 
			 , cocirc_lineage_summary = .lineage_summary( ta )
			 , region_summary = .region_summary( tu )
			 , tips = paste( tu, collapse = '|' )
			 , stringsAsFactors=FALSE
			)
			if ( u %in% report_nodes ) { # print progress 
				i <- which( nodes == u )
				message(paste( 'Progress' , round(100*i / length( nodes )),  '%') ) 
			}
			X
		}, mc.cores = ncpu )
	)
	Y <- Y[ order( Y$logistic_growth_rate , decreasing=TRUE ), ] 
	
	dir.create(output_dir, showWarnings = FALSE)
	ofn1 = glue( paste0( output_dir, '/scanner-{max_date}.rds' ) )
	ofn2 = glue( paste0( output_dir, '/scanner-{max_date}.csv' ) )
	ofn3 = glue( paste0( output_dir, '/scanner-env-{max_date}.rds' ) )
	saveRDS( Y , file=ofn1   )
	write.csv( Y , file=ofn2, quote=FALSE, row.names=FALSE )
	message('saving image ... ' ) 
	e0 = list(  Y = Y 
	  , descendantSids = descendantSids
	  , ancestors = ancestors
	  , sts = sts 
	  , tre = tre
	  , descendantTips = descendantTips
	  , descendants = descendants 
	  , nodedata = nodedata
	  , ndesc = ndesc
	  , descsts = descsts 
	  , clade_age = clade_age
	  , max_desc_time = max_desc_time
	  , min_desc_time = min_desc_time
	  , amd = amd[ , c('sequence_name', 'sample_time', 'sample_date', 'region') ] 
	)  
	saveRDS(e0, file=ofn3)
	message( glue( 'Data written to {ofn1} and {ofn2} and {ofn3}. Returning data frame invisibly.'  ) )
	invisible(Y) 
}

