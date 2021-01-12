#' Compute growth rate statistic for nodes in tree with >x descendants 
#' 
#' Only counts pillar 1 data; comparison sample is matched by time and adm2
#' 
#' @param min_descendants Clade must have at least this many tips 
#' @param max_descendants Clade can have at most this many tips 
#' @param min_date Only include samples after this data 
#' @param max_date Only include samples before and including this date
#' @param ncpu number cpu for multicore ops 
#' @export 
scanner_logistic_growth_rate <- function(min_descendants = 50 , max_descendants = 20e3, min_date = NULL, max_date = NULL , ncpu = 8)
{
	library( treeio ) 
	library( ape ) 
	library( variantAnalysis )
	library( lubridate )
	
	max_time <- Inf 
	if ( !is.null( max_date )){
		max_time <- decimal_date( max_date )
	} else{
		max_date <- Sys.Date()
	}
	
	min_time <- -Inf 
	if (!is.null( min_date ))
		min_time <- decimal_date( min_date )
	
	# play data 
	tre = read.tree(list.files(  '../phylolatest/trees' , patt = 'cog_global_.*_tree.newick', full.names=TRUE) )
	# load coguk algn md 
	amd <- read.csv( list.files(  '../phylolatest/alignments/' , patt = 'cog_[0-9\\-]+_metadata.csv', full.names=TRUE) 
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
	
	.get_comparator_ancestor <- function(u){
		nu = ndesc[u] 
		asu = ancestors[[u]]
		for ( a in asu ){
			na = ndesc[a]
			nc = na - nu 
			if ( nc >= nu )
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
			X = data.frame( node_number = u
			 , parent_number = ifelse( is.null(a), NA, a )
			 , most_recent_tip = as.Date( date_decimal( max( sts[ itu ]  ) ) )
			 , least_recent_tip = as.Date( date_decimal( min( sts[itu]  ) ) )
			 , cluster_size = length( tu )
			 , logistic_growth_rate = .logistic_growth_stat ( u ) 
			 , lineage = paste( unique(amd$lineage[ match( tu, amd$sequence_name)]) , collapse = '|' )
			 , tips = paste( tu, collapse = '|' )
			 , stringsAsFactors=FALSE
			)
			
			cat( paste( Sys.time() , u, X$cluster_size, X$logistic_growth_rate , '\n' ))
			X
		}, mc.cores = ncpu )
	)

	saveRDS( Y , file=glue( 'logisticGrowthStat-{max_date}.rds' )  )
	write.csv( Y , file=glue( 'logisticGrowthStat-{max_date}.csv' ) , quote=FALSE, row.names=FALSE )
	invisible(Y) 
}
