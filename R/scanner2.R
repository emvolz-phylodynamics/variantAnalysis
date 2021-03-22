#' Compute scanner statistics for nodes in tree. 
#'
#' Takes standard inputs in the form of a rooted phylogeny and data frame with required metadata (see below). 
#' Computes a logistic growth rate and a statistic for outlying values in a molecular clock (root to tip regression).
#' Comparison sample is matched by time and region.
#' 
#' TODO: include patientid in metadata & deduplicate 
#' add optional 'sample_weight' to metadata, use in regression models 
#' 
#' @param tre A rooted phylogeny in ape::phylo or treeio::treedata form. The latter may include internal node labels and be read by '.read.beast'. 
#' @param amd A data frame containing required metadata for each tip in tree: sequence_name, sample_date, region. Optional metadata includes: pillar_2, sample_time(numeric), country. 
#' @param min_descendants Clade must have at least this many tips 
#' @param max_descendants Clade can have at most this many tips 
#' @param min_date Only include samples after this data 
#' @param max_date Only include samples before and including this date
#' @param ncpu number cpu for multicore ops 
#' @param output_dir Path to directory where results will be saved 
#' @param include_pillar1 if TRUE (default FALSE), will include Pillar 1 samples when computing stats
#' @param country Optional, only include data from this country 
#' @param num_ancestor_comparison When finding comparison sample for molecular clock stat, make sure sister clade has at least this many tips 
#' @param factor_geo_comparison  When finding comparison sample based on geography, make sure sample has this factor times the number within clade of interest
#' @param generation_time_scale Approximate generation time in years 
#' @export 
scanner2 <- function(tre
 , amd
 , min_descendants = 30 , max_descendants = 20e3, min_date = NULL, max_date = NULL , ncpu = 8
 , output_dir = '.' 
 , include_pillar1 = FALSE 
 , country = NULL 
 , num_ancestor_comparison = 500
 , factor_geo_comparison = 5
 , Tg = 6.5/365
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
		amd$sts <- amd$sample_time 
	}
	if ( !('country' %in% colnames(amd)) )
		amd$country = 'country_not_specified'
	if ( !('pillar_2' %in% colnames(amd)) ){
		amd$pillar_2 = 'True'
	}
	# only data within country 
	if ( !is.null( country ))
		amd <- amd[ amd$country == country, ]
	# exclude p1 
	if ( !include_pillar1 )
	{
		amd <- amd[ amd$pillar_2 == 'True' , ]
	}
	
	# sample time 
	amd <- amd [ (amd$sts >= min_time) & (amd$sts <= max_time) , ] 
	sts <- setNames( amd$sts , amd$sequence_name )
	
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
		# descendants; note also counts mrca node 
		descendants <- lapply( 1:(n+nnode), function(u) u )
		for (ie in 1:nrow(poedges)){
			a <- poedges[ie,1]
			u <- poedges[ie,2]
			descendants[[a]] <- c( descendants[[a]], descendants[[u]] )
		}
		
		# descendant tips
		descendantTips <- lapply( 1:(n+nnode), function(u) c(NA) )
		## tip labels, only including names in amd which meet inclusion criteria 
		descendantSids <- lapply( 1:(n+nnode), function(u) c(NA) )
		for (u in 1:n){
			#descendantTips[[u]] <- u 
			descendantSids[[u]] <- tre$tip.label[u] 
		}
		for (ie in 1:nrow(poedges)){
			a <- poedges[ie,1]
			u <- poedges[ie,2]
			descendantSids[[a]] <- c( descendantSids[[a]], descendantSids[[u]] )
		}
		## remove NA sids (not in amd)
		for (u in (n+1):(n+nnode)) {
			descendantSids[[u]] <- na.omit( descendantSids[[u]] )
			if ( length( descendantSids[[u]] ) > 0 ){
				descendantTips[[u]] <- match(descendantSids[[u]], tre$tip.label ) 
			}
		}
		
		# ancestors
		ancestors <- lapply( 1:(n+nnode), function(u) c() )
		for (ie in 1:nrow(preedges)){
			a <- preedges[ie,1]
			u <- preedges[ie,2]
			ancestors[[u]] <- c( ancestors[[a]], a)
		}
		
		# only counting tips in amd 
		## number tips desc 
		ndesc <- sapply( 1:(n+nnode), function(u) length( descendantSids[[u]] ) )
		## vector samp time desc 
		descsts = lapply( 1:(n+nnode), function(u) sts[ na.omit( descendantSids[[u]] )  ]  )
	}
	print(paste('Derived lookup variables', Sys.time()) ) 
	
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
			#browser() 
			message('Failed to find comparator ancestor with more tips. Returning NA. Node: ', u) 
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
		 stu = descsts [[ u ]]
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
	
	.logistic_growth_stat <- function(u, generation_time_scale = Tg)
	{
		ta = .get_comparator_sample(u) 
		if ( is.null( ta ))
			return(0)
		tu = descendantSids[[u]] 
		sta = sts[ ta ]
		stu = descsts [[ u ]]
		X = data.frame( time = c( sta, stu ), type = c( rep('control', length(ta)), rep('clade',length(tu)) ) )
		X = na.omit( X ) 
		m = glm( type=='clade' ~ time, data = X, family = binomial( link = 'logit' ))
		s = summary( m ) 
		rv = unname( coef( m )[2] * generation_time_scale ) 
		p = NA 
		if ( is.na( rv )){
			message( 'NA growth stat, node: ', u  )	
			#browser()
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
		ta = tryCatch( 
			setdiff( descendantSids[[a]] , tu )
		, error = function(e) browser() )
		sta = sts[ ta ]
		stu = descsts [[ u ]]
		
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
	
	# main 
	## collect nodes with more than x descedants 
	nodes = which(  (ndesc >= min_descendants)   &   (ndesc <= max_descendants) )
	Y = do.call( rbind, 
		parallel::mclapply( nodes , function(u){
			tu = descendantSids[[u]]
			lgs = .logistic_growth_stat ( u )
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
			 , lineage = paste( unique(amd$lineage[ match( tu, amd$sequence_name)]) , collapse = '|' )
			 , tips = paste( tu, collapse = '|' )
			 , stringsAsFactors=FALSE
			)
			X
		}, mc.cores = ncpu )
	)
	
	dir.create(output_dir, showWarnings = FALSE)
	ofn1 = glue( paste0( output_dir, '/scanner-{max_date}.rds' ) )
	ofn2 = glue( paste0( output_dir, '/scanner-{max_date}.csv' ) )
	ofn3 = glue( paste0( output_dir, '/scanner-env-{max_date}.rds' ) )
	saveRDS( Y , file=ofn1   )
	write.csv( Y , file=ofn2, quote=FALSE, row.names=FALSE )
	message('saving image ... ' ) 
	e0 = list( descendantSids = descendantSids, ancestors = ancestors, sts = sts , tre = tre, descendantTips = descendantTips, descendants = descendants , Y = Y 
	  , nodedata = nodedata
	)  
	saveRDS(e0, file=ofn3)
	message( glue( 'Data written to {ofn1} and {ofn2} and {ofn3}. Returning data frame invisibly.'  ) )
	invisible(Y) 
}

