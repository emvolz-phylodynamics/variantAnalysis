#' Compute scanner statistics for nodes in tree. 
#'
#' Takes standard inputs in the form of a rooted phylogeny and data frame with required metadata (see below). 
#' Computes a logistic growth rate and a statistic for outlying values in a molecular clock (root to tip regression).
#' Comparison sample is matched by time and region.
#' 
#' add optional 'sample_weight' to metadata, use in regression models 
#' 
#' @param tre A rooted phylogeny in ape::phylo or treeio::treedata form. The latter may include internal node labels and be read by '.read.beast'. 
#' @param amd A data frame containing required metadata for each tip in tree: sequence_name, sample_date, region. Optional metadata includes: pillar_2, sample_time(numeric), country,  lineage. 
#' @param min_descendants Clade must have at least this many tips 
#' @param max_descendants Clade can have at most this many tips 
#' @param min_date Only include samples after this data 
#' @param max_date Only include samples before and including this date
#' @param min_cluster_age_yrs Only include clades that have sample tips that span at least this value
#' @param min_blen Only compute statistics for nodes descended from branches of at least this length 
#' @param ncpu number cpu for multicore ops 
#' @param output_dir Path to directory where results will be saved 
#' @param include_pillar1 if TRUE (default TRUE), will include Pillar 1 samples when computing stats
#' @param country Optional, only include data from this country 
#' @param num_ancestor_comparison When finding comparison sample for molecular clock stat, make sure sister clade has at least this many tips 
#' @param factor_geo_comparison  When finding comparison sample based on geography, make sure sample has this factor times the number within clade of interest
#' @param Tg Approximate generation time in years 
#' @param report_freq Print progress for every n'th node 
#' @param compute_treestructure If TRUE, will compute treestructure clusters and ensure other tree stats are computed for these clusters
#' @export
scanner2 <- function(tre
 , amd
 , min_descendants = 30 
 , max_descendants = 20e3
 , min_cluster_age_yrs = 1.5/12
 , min_date = NULL
 , max_date = NULL 
 , min_blen = min( tre$edge.length[ tre$edge.length > 0 ] )/2
 , ncpu = 8
 , output_dir = '.' 
 , include_pillar1 = TRUE 
 , country = NULL 
 , num_ancestor_comparison = 500
 , factor_geo_comparison = 5
 , Tg = 6.5/365
 , report_freq = 50 
 , compute_treestructure = FALSE
 , lsd_args = list(	givenRate = .00075
		, seqLen = 29e3
		, variance = 0 
		, verbose = TRUE 
	) 
  , ts_args = list( minCladeSize= 50
	, level = .001
  )
  , root_on_tip = 'Wuhan/WH04/2020'
  , root_on_tip_sample_time = 2020 
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
		if ( compute_treestructure ){
			stop( 'Can not compute treestructure statistics when input tree is not an ape::phylo' )
		}
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
	stopifnot( all(tre$tip.label %in% amd$sequence_name) )
	if ( !('sample_time' %in% colnames(amd))){
		amd$sample_time = decimal_date (amd$sample_date)
	}
	amd$sts <- amd$sample_time 
	if (!('lineage' %in% colnames( amd ))){
		amd$lineage <- 'lineage_not_provided'
	}
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
	
	
	# retain only required variables 
	amd = amd[ , c('sequence_name', 'sample_time', 'sample_date', 'region', 'lineage') ] 
	amd <- amd [ amd$sequence_name %in% tre$tip.label ,  ] 
	
	# prune tree
	if ( !( root_on_tip %in% amd$sequence_name)){
		amd <- rbind( amd, data.frame( 
			sequence_name = root_on_tip
			, sample_time = root_on_tip_sample_time
			, sample_date = as.Date( lubridate::date_decimal( root_on_tip_sample_time ) )
			, region = NA
			, lineage = NA
		))
	}
	tre <- keep.tip( tre, intersect( tre$tip.label , amd$sequence_name )  )
	if (  !( root_on_tip %in% tre$tip.label ) )  {
		stop('Outgroup sequence missing from input tree.')
	}
	tr2 = root( tre, outgroup= root_on_tip, resolve.root = TRUE )
	
	sts <- setNames( amd$sample_time[ match( tr2$tip.label, amd$sequence_name ) ] , tr2$tip.label )

	
	## treat tip labs not in amd differently 
	if ( compute_treestructure ) {
		# prune the tree to match amd and make rooted
		# date tree 
		# make bifurcating 
		# run treestructure 
		library( treestructure ) 
		library (Rlsd2 )
		
		message( paste( Sys.time() , 'Starting time tree' ))
		lsd_args$inputTree = tr2 
		lsd_args$inputDate = sts 
		lt1 = do.call( Rlsd2::lsd2, lsd_args ) 
		message( paste( Sys.time() , 'Made time tree' )) #45 mins 
		
		otr = as.phylo( lt1$dateNexusTreeFile ) # not a date or a nexus or a file, but it is a tree
		otr1 = root( otr, root_on_tip , resolve.root=TRUE) # reroot tree and extract relevant clade 
		iog = which(otr1$tip.label==root_on_tip)
		rootnode = otr1$edge[ otr1$edge[,2]==iog ,1 ]
		uv = otr1$edge[  otr1$edge[,1] == rootnode , 2]
		u <- setdiff( uv, iog )
		# overwriting 'tre'
		tre = multi2di( extract.clade( otr1, u, root.edge = 1 )) # the timed tree excluding outgroup, made bifurcating 
		message( paste( Sys.time() , 'Made bifurcating time tree ' ))
		ts_args$ncpu = ncpu 
		ts_args$tre = tre 
		ts <- do.call( treestructure::trestruct.fast, ts_args )
		message( paste( Sys.time() , 'Computed treestructure clusters' ))
	} 
	
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
		## time dep growth 
		m1 = mgcv::gam( type=='clade' ~ s(sample_time, bs = "bs", k = 4, m=1) , family = binomial(link='logit') , data = X)
		tout = seq( min(sample_time) , max(sample_time), length=5)
		tout1 = tout[4] + diff(tout)[1]/2 
		dlo = diff( predict( m1, newdata = data.frame( sample_time = c( tout1 , max(tout))) ) 
		r=dlo*7 / ((max(tout) - tout1)*365) 
		c( lgr = rv, lgrp = p, gam_r =r , dAIC = AIC(m1) - AIC(m) )
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
	
	# Compute 'proportionality' statistics for node (u) and given variable (var)
	# e.g. if sample is from vaccine breakthrough, is there higher odds that sample is in clade?
	.var_proportionality_stat <- function (u, var = 'breakthrough2'
	 , f = clade ~ time + var 
	 , form_index = 3 
	 , value=TRUE, ta = NULL) 
	{
		if (is.null(ta)) 
			ta = .get_comparator_sample(u)
		if (is.null(ta)) 
			return(0)
		tu = descendantSids[[u]]
		sta = sts[ta]
		stu = sts[tu]
		X = data.frame(tip = c( ta, tu )
					, time = c(sta, stu)
					, type = c(rep("control", length(ta)), rep("clade", length(tu)))
		)
		X$var <- (amd[[var]][match(X$tip, amd$sequence_name)]==value)
		X$clade <- (X$type == 'clade' )
		
		#m = glm(  var ~ time + (type=='clade'), data = X, family = binomial(link = "logit"))
		m = glm(  f, data = X, family = binomial(link = "logit"))
		s = summary(m)
		rv = unname(coef(m)[form_index])
		p = NA
		if (is.na(rv)) {
			message("NA proportionality stat, node: ", u)
		}
		else {
			p = s$coefficients[form_index, 4]
		}
		c(rv, p)
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
	nodes = which(  (ndesc >= min_descendants)   &   (ndesc <= max_descendants) & (clade_age >= min_cluster_age_yrs)   )
	nodes_blenConstraintSatisfied <- tre$edge[ tre$edge.length > min_blen , 2]
	nodes <- intersect( nodes, nodes_blenConstraintSatisfied )
	if (compute_treestructure) {
		nodes <- union( ts$cluster_mrca, nodes )
	}
	report_nodes <- nodes[ seq(1, length(nodes), by = report_freq) ] #progress reporting 
	node_ancestors <- do.call( c, lapply( nodes, function(u) ancestors[[u]] ))
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
			 , lineage_summary = tryCatch( .lineage_summary( tu ) , error = function(e) as.character(e) )
			 , cocirc_lineage_summary = tryCatch( .lineage_summary( ta ), error = function(e) as.character(e))
			 , region_summary = tryCatch( .region_summary( tu ), error = function(e) as.character(e))
			 , external_cluster = !(u %in% node_ancestors ) 
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
	if ( ! all( nodes %in% Y$node_number)){
		stop('Statistics not computed for all nodes. Possible that memory exceeded with ncpu > 1. Try with ncpu = 1.') 
	}
	Y <- Y[ order( Y$logistic_growth_rate , decreasing=TRUE ), ] 
	
	## add in trestruct stat 
	if ( compute_treestructure ){
		zdf = data.frame( node_number = ts$cluster_mrca , treestructure_z = ts$cluster_z)
		Y = merge( Y, zdf, all.x = TRUE, sort = FALSE  , by = 'node_number' )
	} else{
		Y$treestructure_z = NA 
	}
	
	dir.create(output_dir, showWarnings = FALSE)
	ofn1 = glue( paste0( output_dir, '/scanner-{max_date}.rds' ) )
	#ofn2 = glue( paste0( output_dir, '/scanner-{max_date}.csv' ) )
	ofn3 = glue( paste0( output_dir, '/scanner-env-{max_date}.rds' ) )
	saveRDS( Y , file=ofn1   )
	#write.csv( Y , file=ofn2, quote=FALSE, row.names=FALSE )
	message('saving image ... ' ) 
	
	# save internal variables and functions 
	e0 = environment() 
	
	saveRDS(e0, file=ofn3)
	message( glue( 'Data written to {ofn1} and {ofn3}. Returning data frame invisibly.'  ) )
	invisible(Y) 
}



#'
#'
#' @export 
get_clusternode_gam <- function( u0, Y, e0, ... ){
	tu = strsplit( Y$tips[ Y$node_number == u0 ], '\\|' )[[1]]
	ta = e0$.get_comparator_sample( u0, e0, nX = 5 ) 
	s <- e0$amd 
	s <- s[ s$sequence_name  %in% c(tu, ta ), ]
	s$genotype <- 'non-Cluster'
	s$genotype[s$sequence_name %in% tu] <- paste( 'Cluster', u0 )
	library( ggplot2 )
	vfd <- variantAnalysis::variable_frequency_day(s, variable = "genotype", value = paste( 'Cluster', u0 ), ...) 
	vfd
}


#'
#'
#' @export 
get_clusternode_gam_maps <- function(u0, Y, e0 )
{
	s1 <- e0$amd
	if (  ( 'epi_week' %in% colnames(s1))  &  (!('epiweek' %in% colnames(s1) ) ) ){
		s1$epiweek<- s1$epi_week 
	}else if  (  !( 'epi_week' %in% colnames(s1))  &  (!('epiweek' %in% colnames(s1) ) ) ){
		ew <- lubridate::epiweek( as.Date(s1$sample_date)  )
		i <- with( s1 ,  (year(as.Date(s1$sample_date)) > 2020)  &  (ew < 53)  )
		ew[ i ] <- 53 + ew[i] 
		s1$epiweek <- ew 
	}
	tu = strsplit( Y$tips[ Y$node_number == u0 ], '\\|' )[[1]]
	
	library(raster)
	library( spdep )
	library( cowplot ) 
	library( ggplot2 )

	shp <- shapefile(system.file('Local_Authority_Districts__December_2019__Boundaries_UK_BUC.shp', package='variantAnalysis') )
	shp <- shp[ shp$lad19cd %in% unique( s1$ltlacd ), ]
	shpdf <- droplevels(as(shp, 'data.frame'))
	nb <- poly2nb(shp, row.names = shpdf$lad19cd)
	names(nb) <- attr(nb, "region.id")

	s1$VUI <- FALSE
	s1$VUI[s1$sequence_name %in% tu] <- TRUE
	s1$n <- 1 
	s2 = merge( aggregate( VUI ~ epiweek*ltlacd,  data = s1 , sum )
	 , aggregate( n ~ epiweek*ltlacd,  data = s1 , length ) 
	 , by = c('epiweek' , 'ltlacd' )
	)
	toadd <- setdiff( as.factor( s$ltlacd ) , s2$fac_ltlacd )
	s2 <- rbind( s2
		, data.frame( epiweek = max(s2$epiweek), ltlacd = toadd, VUI=0, n = 0 ) 
	)
	s2$fac_ltlacd = as.factor( s2$ltlacd )
	s2 <- na.omit( s2 )
	m1 = mgcv::gam( cbind(VUI, n - VUI) ~ s(epiweek) + s(fac_ltlacd, bs='mrf', xt=list(nb=nb), k = 50), family = binomial(link='logit'), data = s2,  method = 'REML', ctrl=gam.control(nthreads = 6))

	s2$logodds = predict( m1 )
	s2s <- split( s2, s2$ltlacd )
	 
	shp0 <- shapefile(system.file('Local_Authority_Districts__December_2019__Boundaries_UK_BUC.shp', package='variantAnalysis') )
	shp0 <- fortify(shp0, region = 'lad19cd')
	maxfreq <- max( with( s2,  exp( logodds ) / ( 1 + exp( logodds ) )  ))
	plmap = function(EPIWEEK, cscale = FALSE ){
		shp = shp0 
		idf2 <- do.call( rbind, lapply( s2s, function(x)  x[x$epiweek == EPIWEEK  , ] ))
		idf2$id <- idf2$ltlacd
		shp <- merge(shp, idf2, by = 'id', all.x = TRUE)
		shp <- dplyr::arrange(shp, order)
		shp$frequency_VUI <- exp( shp$logodds ) / ( 1 + exp( shp$logodds ) )
		shp1 = shp 
		p0 <- ggplot(data = shp1, aes(x = long, y = lat, group = group, fill = frequency_VUI)) + geom_polygon() + coord_equal() + theme_void() + theme(legend.title = element_blank()) + ggtitle( glue( 'Epi week {EPIWEEK}' )) + scale_fill_gradient(low ='lightgrey', high ='red' ,limits = c(0,maxfreq))
		
		# '#f8766dff' , '#00bfc4ff'
		#  '#cc993373' , '#3366cc64' 
		if (!cscale) 
			p0 <- p0 + theme( legend.pos = 'none' )
		p0
	}

	MAXWEEK = max( s1$epiweek[s1$VUI] )-1
	MINWEEK = MAXWEEK - 8
	pls = lapply( MINWEEK:MAXWEEK, plmap )
	pls[[6]] <- plmap( (MINWEEK:MAXWEEK)[6], TRUE ) 

	p = plot_grid( plotlist = pls , nrow = 3 )
	lwp = plmap( MAXWEEK, TRUE ) + ggtitle(paste('Cluster', u0)) 
	lwp2 = lwp + theme( legend.key.size = unit(.25, "in"), legend.key.width = unit(0.1,"in") )

	list( lwp2, p )
}



#'
#'
#' @export 
get_clusternode_mlesky <- function( u=406318 , scanner_env=readRDS("scanner-env-2021-03-03.rds")) 
{
  library(mlesky)
  library( treedater )
  
  e1 = as.environment( scanner_env )
  attach( e1 )
  
  mr = 5.9158E-4
  utre <- keep.tip(tre, descendantSids[[u]])
  sample_times <- sts[utre$tip.label]
  
  tr <- di2multi(utre, tol = 1e-05)
  tr = unroot(multi2di(tr))
  tr$edge.length <- pmax(1/29000/5, tr$edge.length)
  tr3 <- dater(unroot(tr), sts[tr$tip.label], s = 29000, omega0 = mr, 
               numStartConditions = 0, meanRateLimits = c(mr, mr + 1e-6), ncpu = 6)
  
  msg <- mlskygrid(tr3, tau = NULL, tau_lower=.001, tau_upper = 10 , sampleTimes = sts[tr3$tip.label] , 
                   res = 10, ncpu = 3)
  
  list( mlesky = msg, timetree = tr3, tree = tr  )
}



#' Computes mutations for sequences within specified clusters
#'
#' Also derives defining mutations by contrasting with closest ancestor
#'
#' @param scanner_env output of scanner
#' @param nodes integer vector of nodes for which results will be computed. Computes for all nodes in Y if omitted
#' @param mutdata path to file or data frame with variant information
#' @param mut_variable Name of variable in data frame (spec by mutfn) which contains genetic variant data
#' @param ncpu number of CPUs for parallel processing 
#' @export 
cluster_muts = function( 
  scanner_env 
 , nodes = NULL
 , mutdata = list.files( patt = 'cog_global_[0-9\\-]+_mutations.csv' , full.name=TRUE , path = '../phylolatest/metadata/' )
 , mut_variable = c( 'variants', 'mutations' )
 , min_seq_contrast = 1 # sequences in ancestor clade
 , overlap_threshold = .9
 , ncpu = 1
 )
{
	
	e1 = as.environment( scanner_env )
	attach( e1 )
	
	mut_variable = mut_variable[1] 
	if ( !('tips' %in% colnames( Y)))
		stop('Input data is missing column `tips`')
	
	
	if (is.null( nodes ))
		nodes <- Y$node_number
	sdnodes = setdiff( nodes, Y$node_number )
	if ( length( sdnodes ) > 0 ){
		warning( paste('Some input nodes not found in scanner table', sdnodes, collapse = ', ') ) 
	}
	nodes1 <- setdiff( nodes, sdnodes ) #intersect( Y$node_number, nodes )
	#Ygr1 = Y[ Y$node_number %in% nodes , ] 
	Ygr1 <- Y[ match( nodes1, Y$node_number ), ]
	
	if ( is.character( mutdata )){
		mdf = read.csv( mutdata, stringsAs=FALSE )
	} else if (is.data.frame(mutdata)) {
		mdf = mutdata
	}else if ( !is.data.frame( mutdata )  ) {
		stop('*mutdata* must be a data frame or a path to a csv ')
	}
	
	cms = parallel::mclapply( 1:nrow(Ygr1) , function(ku) {
		tipsku = strsplit( Ygr1$tips[ku], split = '\\|')[[1]] 
		mdf.u = mdf[ mdf$sequence_name %in% tipsku, ]
		u = nodes[ku] 
		
		## find comparator ancestor 
		for ( a in rev( ancestors[[u]] ) ){ ## NOTE ancestors goes in order of (tree root) -> u , so rev to go from u to tree root 
			sdt = setdiff( descendantTips[[a]] , descendantTips[[u]] ) 
			if ( length( sdt ) >= min_seq_contrast )
				break 
		}
		asids = setdiff( descendantSids[[a]] ,  descendantSids[[u]] )
		mdf.a = mdf[ mdf$sequence_name %in% asids   , ]
		if ( nrow( mdf.a ) == 0  |  nrow( mdf.u ) == 0 )
			return( list(defining = NA, all = NA ) )
		vtabu = sort( table( do.call( c, strsplit( mdf.u[[mut_variable]], split='\\|' )  ) ) / nrow( mdf.u ) )
		vtaba = sort( table( do.call( c, strsplit( mdf.a[[mut_variable]], split='\\|' )  ) ) / nrow( mdf.a ) )
		
		umuts = names( vtabu[ vtabu > overlap_threshold ] )
		defining_muts = setdiff( names( vtabu[ vtabu > overlap_threshold ] )
		 , names(vtaba[ vtaba > overlap_threshold ]) )
		list(defining=defining_muts, all=umuts  ) 
	}, mc.cores = ncpu)
	
	detach( e1 )
	
	names( cms ) <- Ygr1$node_number 
	cms
}


#' 
#'
#' @export 
condense_clusters <- function(  scanner_env, threshold_growth = .5 , candidate_nodes = NULL){
  e1 = as.environment( scanner_env )
  attach( e1 )
  if ( is.null( candidate_nodes )){
    candidate_nodes = cnodes = Y$node_number[ Y$logistic_growth_rate >= threshold_growth ]
  }
  stopifnot( length( cnodes ) > 0 )
  keep <- sapply( cnodes, function(u){
    tu = na.omit( descendantTips[[u]] )
    x = sapply( setdiff(cnodes,u) , function(a) {
      ta = na.omit(  descendantTips[[a]] )
      alltuta = (all(tu %in% ta))
      y = (length(ta) > length(tu))  &  alltuta
      z = (length(ta)==length(tu)) & alltuta 
      if (z & (a %in% ancestors[[u]])) {
        return(FALSE) 
      } else if (z & !(a %in% ancestors[[u]]) ){
        return(TRUE)
      }
      y
    })
    all( !x )
  })
  keepnodes <- cnodes[ keep ]
  detach( e1 ) 
  keepnodes 
}


#~ treestructure_scan <- function( e0, omega = .00075, s = 29e3, minCladeSize = 50, level = .001,  ncpu = 6, ...) 
#~ {
#~ 	library( treestructure ) 
#~ 	library( Rlsd2 )
#~ 	library( ape ) 
	
#~ 	attach( e0 )
#~ 	message( paste( Sys.time() , 'Starting' ))
#~ 	tr1 <- multi2di( tre )
#~ 	message( paste( Sys.time() , 'Made bifurcating tree ' ))
#~ 	sts <- setNames( amd$sample_time[ match( tr1$tip.label, amd$sequence_name ) ] , tr1$tip.label )
#~ 	if ( any ( is.na( sts ))){
#~ 		i = which( is.na( sts ))
#~ 		tr1 <- drop.tip( tr1, names(sts)[i] )
#~ 		sts <- sts[ tr1$tip.label ]
#~ 	} 
#~ 	lt1 = lsd2( inputTree=tr1, inputDate= sts, givenRate = omega, seqLen = s, constraint=FALSE) 
#~ 	message( paste( Sys.time() , 'Made time tree' ))
#~ 	otr = as.phylo( lt1$dateNexusTreeFile ) 
#~ 	ts = trestruct.fast( otr, minCladeSize = minCladeSize , level = level,  ncpu = ncpu, ...) 
#~ 	message( paste( Sys.time() , 'Computed clusters' ))
#~ 	detach( e0 )
	
#~ 	ts
#~ }
