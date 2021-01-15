
#' Compute scannter statistics for nodes in tree with >x descendants 
#' 
#' Computes a logistic growth rate and a statistic for outlying values in a molecular clock (root to tip regression).
#' Only counts pillar 1 data; comparison sample is matched by time and adm2.
#' 
#' @param treenexfn File name of tree with metadata in node labels
#' @param min_descendants Clade must have at least this many tips 
#' @param max_descendants Clade can have at most this many tips 
#' @param min_date Only include samples after this data 
#' @param max_date Only include samples before and including this date
#' @param ncpu number cpu for multicore ops 
#' @param path_to_data Path to data with COG tree and COG tree metadata 
#' @export 
scanner <- function(treenexfn = NULL, min_descendants = 30 , max_descendants = 20e3, min_date = NULL, max_date = NULL , ncpu = 8
 , path_to_data = '/cephfs/covid/bham/results/phylogenetics/latest/'
 )
{
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
	if ( is.null( treenexfn ) ){
		tre = read.tree(list.files(  paste0(path_to_data, '/trees') , patt = 'cog_global_.*_tree.newick', full.names=TRUE) )
		nodedata = NULL
	} else{
		library( treeio ) 
		trd = read.beast( treenexfn )
		tre = get.tree( trd )
		nodedata = as.data.frame( get.data( trd ) )
		rownames( nodedata ) <- nodedata$node
	}
	
	# load coguk algn md 
	amd <- read.csv( list.files(  paste0( path_to_data , '/alignments/') , patt = 'cog_[0-9\\-]+_metadata.csv', full.names=TRUE) 
	  , stringsAs=FALSE )
	amd$sample_time = decimal_date ( as.Date( amd$sample_date ))
	
	# exclude global 
	tre = keep.tip( tre, tre$tip.label[ tre$tip %in% amd$sequence_name[ amd$country=='UK' ] ]) 
	# exclude p1 
	lhls <- paste0( '.*/', c( 'MILK', 'ALDP', 'QEUH', 'CAMC'), '.*')
	lhpatt = paste( lhls, collapse = '|' )
	tre = keep.tip(tre,  tre$tip.label[ grepl( tre$tip.label, patt = lhpatt ) ] ) 
	amd <- amd [ amd$sequence_name %in% tre$tip.label , ]
	
	# sample time 
	sts <- decimal_date ( as.Date( amd$sample_date[ match( tre$tip.label, amd$sequence_name ) ] ) ) 
	names(sts) <- tre$tip.label 
	todrop <- names(sts) [ sts < min_time | sts > max_time ]
	tre <- drop.tip( tre, todrop )
	sts <- sts [ tre$tip.label ]
	amd <- amd [ match( tre$tip.label, amd$sequence_name ) , ]
	amd <- amd[ !is.na( amd$sequence_name ) , ]
	
	# root to tip 
	ndel <- node.depth.edgelength( tre ) 
	
	# data structures to quickly look up tree data 
	# copied/adapted from treestructure 
	{
		n <- Ntip( tre ) 
		nnode = Nnode( tre )
		poedges <- tre$edge[ ape::postorder( tre ), ]
		preedges <- tre$edge[ rev( ape::postorder( tre )), ]
		
		# PRECOMPUTE for each node 
		# descendants; note also counts mrca node 
		descendants <- lapply( 1:(n+nnode), function(u) u )
		for (ie in 1:nrow(poedges)){
			a <- poedges[ie,1]
			u <- poedges[ie,2]
			descendants[[a]] <- c( descendants[[a]], descendants[[u]] )
		}
		
		# descendant tips
		descendantTips <- lapply( 1:(n+nnode), function(u) c() )
		for (u in 1:n)
		  descendantTips[[u]] <- u 
		for (ie in 1:nrow(poedges)){
			a <- poedges[ie,1]
			u <- poedges[ie,2]
			descendantTips[[a]] <- c( descendantTips[[a]], descendantTips[[u]] )
		}
		
		# descendant internal
		descInternal <- lapply( 1:(n+nnode), function(u) c() )
		for (u in (n+1):(n+nnode))
		  descInternal[[u]] <- u 
		for (ie in 1:nrow(poedges)){
			a <- poedges[ie,1]
			u <- poedges[ie,2]
			descInternal[[a]] <- c( descInternal[[a]], descInternal[[u]] )
		}
		
		# ancestors
		ancestors <- lapply( 1:(n+nnode), function(u) c() )
		for (ie in 1:nrow(preedges)){
			a <- preedges[ie,1]
			u <- preedges[ie,2]
			ancestors[[u]] <- c( ancestors[[a]], a)
		}
		
		# number tips desc 
		ndesc <- sapply( 1:(n+nnode), function(u) length( descendantTips[[u]] ) )
		
		# vector samp time desc 
		descsts = lapply( 1:(n+nnode), function(u) sts[ descendantTips[[u]] ]  )
	}
	
	.get_comparator_ancestor <- function(u, num_comparison = 500)
	{
		nu = ndesc[u] 
		asu = ancestors[[u]]
		for ( a in asu ){
			na = ndesc[a]
			nc = na - nu 
			if ( nc >= num_comparison )
			{
				break 
			}
		}
		if ( na < nu )
			return (NA) 
		a
	}
	#~ a = .get_comparator_ancestor ( u )
	
	# matched by time and in proportion to adm2 prevalence
	.get_comparator_sample <- function( u , nX = 5 ) 
	{
		# weight for adm2 
		 w = table ( amd$adm2[ match( tre$tip[ descendantTips[[u]]  ] , amd$sequence_name ) ]  )  
		 w = w [ names(w)!='' ]
		 w = w / sum( w ) 
		 
		 nu <- ndesc[u] 
		 tu = tre$tip.label [ descendantTips[[u]] ]
		 stu = descsts [[ u ]]
		 minstu = min(na.omit(stu ))
		 maxstu = max(na.omit(stu) )
		 
		 
		 amd1 = amd[ (amd$sample_time >= minstu) & (amd$sample_time <= maxstu), c('sequence_name', 'sample_date', 'adm2') ]
		 amd1 <- amd1[ !(amd1$sequence_name %in% tu) , ]
		 amd1 <- amd1[ amd1$adm2 %in% names( w ) , ]
		 na = min ( nrow( amd1 ) , nu * nX )
		 if ( na < nu ) 
			return ( NULL )
		 amd1$w = w[ amd1$adm2 ] 
		 ta = sample( amd1$sequence_name, replace=FALSE, size = na , prob = amd1$w ) 
		 ta
	}

	.logistic_growth_stat <- function(u, generation_time_scale = 6.5/365)
	{
		ta = .get_comparator_sample(u) 
		if ( is.null( ta ))
			return(0)
		tu = tre$tip.label[  descendantTips[[u]]  ]
		sta = sts[ ta ]
		stu = descsts [[ u ]]
		X = data.frame( time = c( sta, stu ), type = c( rep('control', length(ta)), rep('clade',length(tu)) ) )
		m = glm( type=='clade' ~ time, data = X, family = binomial( link = 'logit' ))
		summary( m ) 
		unname( coef( m )[2] * generation_time_scale ) 
	}
	
	# log median p-value of rtt predicted divergence of tips under u
	.clock_outlier_stat <- function(u)
	{
		a = .get_comparator_ancestor(u)
		if ( is.na( a ))
			return( NA ) 
		
		iu = descendantTips[[u]]
		tu = tre$tip.label [ iu ]
		ia = setdiff( descendantTips[[a]], iu )
		ta = tre$tip.label [ ia ]
		sta = sts[ ta ]
		stu = descsts [[ u ]]
		
		ndelu = ndel[ iu ] 
		ndela = ndel[ ia ] 
		m = lm ( ndela ~ sta ) 
		r2 = summary( m )$r.squared
		oosp = predict(m, newdata =  data.frame(sta = unname( stu )) )  -  ndelu
		mean( oosp^2) / mean( (predict(m) - ndela )^2 )
	}
	
	# main 
	## collect nodes with more than x descedants 
	nodes = which(  (ndesc >= min_descendants)   &   (ndesc <= max_descendants) )
	#~ node_number	parent_number	most_recent_tip	least_recent_tip	day_range	persistence	recency	age	tip_count	uk_tip_count	uk_child_count	uk_chain_count	identical_count	divergence_ratio	mean_tip_divergence	stem_length	growth_rate	lineage	uk_lineage	proportion_uk	admin0_count	admin1_count	admin2_count	admin0_mode	admin1_mode	admin2_mode	admin1_entropy	admin2_entropy	tips
	Y = do.call( rbind, 
		parallel::mclapply( nodes , function(u){
			itu = descendantTips[[u]] 
			tu = tre$tip.label[itu]
			a = head(ancestors[[u]] , 1 )
			na = ndesc[a]
			X = data.frame( cluster_id = ifelse(is.null(nodedata), as.character(u) ,  nodedata[ as.character(u), 'cluster_id' ] ) 
			 , most_recent_tip = as.Date( date_decimal( max( sts[ itu ]  ) ) )
			 , least_recent_tip = as.Date( date_decimal( min( sts[itu]  ) ) )
			 , cluster_size = length( tu )
			 , logistic_growth_rate = .logistic_growth_stat ( u ) 
			 , clock_outlier = .clock_outlier_stat(u)
			 , lineage = paste( unique(amd$lineage[ match( tu, amd$sequence_name)]) , collapse = '|' )
			 , tips = paste( tu, collapse = '|' )
			 , stringsAsFactors=FALSE
			)
			
			cat( paste( Sys.time() , u, X$cluster_size, X$logistic_growth_rate, X$clock_outlier , '\n' ))
			X
		}, mc.cores = ncpu )
	)
	
	ofn1 = glue( 'scanner1-{max_date}.rds' )
	ofn2 = glue( 'scanner1-{max_date}.csv' )
	saveRDS( Y , file=ofn1   )
	write.csv( Y , file=ofn2, quote=FALSE, row.names=FALSE )
	cat( glue( 'Data written to {ofn1} and {ofn2}. Returning data frame invisibly.\n \n'  ) )
	invisible(Y) 
}
