
#' Compute scannter statistics for nodes in tree with >x descendants 
#' 
#' Computes a logistic growth rate and a statistic for outlying values in a molecular clock (root to tip regression).
#' Only counts pillar 1 data; comparison sample is matched by time and adm2.
#' 
#' @param treenexfn File name of tree with metadata in node labels; if NULL, will instead use latest tree found in path_to_data
#' @param min_descendants Clade must have at least this many tips 
#' @param max_descendants Clade can have at most this many tips 
#' @param min_date Only include samples after this data 
#' @param max_date Only include samples before and including this date
#' @param ncpu number cpu for multicore ops 
#' @param path_to_data Path to data with COG tree and COG tree metadata 
#' @param output_dir Path to directory where results will be saved 
#' @param include_pillar1 if TRUE (default FALSE), will include Pillar 1 samples when computing stats
#' @export 
scanner <- function(treenexfn = NULL, min_descendants = 30 , max_descendants = 20e3, min_date = NULL, max_date = NULL , ncpu = 8
 , path_to_data = '/cephfs/covid/bham/results/phylogenetics/latest/'
 , output_dir = '.' 
 , include_pillar1 = FALSE
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
	if ( is.null( treenexfn ) ){
		tre = read.tree(list.files(  paste0(path_to_data, '/trees') , patt = 'cog_global_.*_tree.newick', full.names=TRUE) )
		nodedata = NULL
	} else{
		library( treeio ) 
		trd = .read.beast( treenexfn )
		tre = get.tree( trd )
		nodedata = as.data.frame( get.data( trd ) )
		rownames( nodedata ) <- nodedata$node
	}
	
	# load coguk algn md 
	amd <- read.csv( list.files(  paste0( path_to_data , '/alignments/') , patt = 'cog_[0-9\\-]+_metadata.csv', full.names=TRUE) 
	  , stringsAs=FALSE )
	amd$sample_time = decimal_date ( as.Date( amd$sample_date ))
	
	# exclude global 
	amd <- amd[ amd$country == 'UK', ]
	# exclude p1 
	if ( !include_pillar1 )
	{
		lhls <- paste0( '.*/', c( 'MILK', 'ALDP', 'QEUH', 'CAMC'), '.*')
		lhpatt = paste( lhls, collapse = '|' )
		amd <- amd[ grepl(amd$sequence_name, patt = lhpatt ) , ] 
	}
	
	# sample time 
	amd$sts <- decimal_date ( as.Date( amd$sample_date ) )
	amd <- amd [ (amd$sts >= min_time) & (amd$sts <= max_time) , ] 
	sts <- setNames( amd$sts , amd$sequence_name )
	amd <- amd[ !is.na( amd$sequence_name ) , ]
	
	## treat tip labs not in amd differently 
	tre$tip.label [ !(tre$tip.label %in% amd$sequence_name) ] <- NA 
	
	# root to tip 
	ndel <- node.depth.edgelength( tre ) 
	
	print(paste('Loaded tree & filtered by inclusion criteria', Sys.time()) )
	
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
	
	.get_comparator_ancestor <- function(u, num_comparison = 500)
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
	
	# matched by time and in proportion to adm2 prevalence
	.get_comparator_sample <- function( u , nX = 5 ) 
	{
		# weight for adm2 
		 w = table ( amd$adm2[ match( descendantSids[[u]]  , amd$sequence_name ) ]  )  
		 w = w [ names(w)!='' ]
		 w = w / sum( w ) 
		 
		 nu <- ndesc[u] 
		 tu =  descendantSids[[u]] 
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
	#~ node_number	parent_number	most_recent_tip	least_recent_tip	day_range	persistence	recency	age	tip_count	uk_tip_count	uk_child_count	uk_chain_count	identical_count	divergence_ratio	mean_tip_divergence	stem_length	growth_rate	lineage	uk_lineage	proportion_uk	admin0_count	admin1_count	admin2_count	admin0_mode	admin1_mode	admin2_mode	admin1_entropy	admin2_entropy	tips
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
	cat('saving image ... \n' ) 
	e0 = list( descendantSids = descendantSids, ancestors = ancestors, sts = sts , tre = tre, descendantTips = descendantTips, descendants = descendants , Y = Y 
	  , nodedata = nodedata
	)  
	saveRDS(e0, file=ofn3)
	cat( glue( 'Data written to {ofn1} and {ofn2} and {ofn3}. Returning data frame invisibly.\n \n'  ) )
	invisible(Y) 
}

#' wrapper for external calls to inner_condense_clusters
#'
#' @export
condense_clusters <- function( scanner_env, threshold_growth = .5 , candidate_nodes = NULL){
  e1 = as.environment( scanner_env )
  attach( e1 )
  keepnodes <- inner_condense_clusters(Y, threshold_growth , candidate_nodes)
  detach( e1 )
  keepnodes 
}

#' internal function calls only
inner_condense_clusters <- function( Y, threshold_growth = .5 , candidate_nodes = NULL){
  if ( is.null( candidate_nodes )){
    candidate_nodes = cnodes = Y$node_number[ Y$logistic >= threshold_growth ]
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
  keepnodes 
}

#' @export 
cluster_muts = function( Y
 , scanner_env 
 , nodes 
 , mutfn = list.files( patt = 'cog_global_[0-9\\-]+_mutations.csv' , full.name=TRUE , path = '../phylolatest/metadata/' )
 , min_seq_contrast = 1 # sequences in ancestor clade
 , overlap_threshold = .9
 )
{
	e1 = as.environment( scanner_env )
	attach( e1 )
	
	Ygr1 = Y[ Y$node_number %in% nodes , ] 
	
	mdf = read.csv( mutfn, stringsAs=FALSE )
	
	cms = lapply( 1:nrow(Ygr1) , function(ku) {
		mdf.u = mdf[ mdf$sequence_name %in% strsplit( Ygr1$tips[ku], split = '\\|')[[1]] , ]
		u = nodes[ku] 
		
		## find comparator ancestor 
		for ( a in ancestors[[u]][-1] ){
			sdt = setdiff( descendantTips[[a]] , descendantTips[[u]] ) 
			if ( length( sdt ) >= min_seq_contrast )
				break 
		}
		asids = setdiff( descendantSids[[a]] ,  descendantSids[[u]] )
		mdf.a = mdf[ mdf$sequence_name %in% asids   , ]

		vtabu = sort( table( do.call( c, strsplit( mdf.u$variants, split='\\|' )  ) ) / nrow( mdf.u ) )
		vtaba = sort( table( do.call( c, strsplit( mdf.a$variants, split='\\|' )  ) ) / nrow( mdf.a ) )
		
		res = setdiff( names( vtabu[ vtabu > overlap_threshold ] )
		 , names(vtaba[ vtaba > overlap_threshold ]) )
		sort( res ) 
	})
	
	detach( e1 )
	
	ress2 = cms 
	if ( length(cms) > 1 ){
		ress2 = lapply( cms, function( x ) setdiff( x, Reduce( intersect, cms )  ) )
	} 
	ress2
}
#~ cluster_muts( Y, ey, nodes = 449569 )



#' matched by time and in proportion to adm2 prevalence
#'
#' @export 
get_comparator_sample <- function( u , scanner_env, nX = 5 ) 
{
	e1 = as.environment( scanner_env )
	attach( e1 )
	 ta <- inner_get_comparator_sample(u, nX)
	 detach( e1 )
	 ta
}

#' internal function calls only, assuming amd, ndesc, descsts are already in the environment
#'
inner_get_comparator_sample <- function( u, nX = 5 ) 
{
  # weight for adm2 
  w = table ( amd$adm2[ match( descendantSids[[u]]  , amd$sequence_name ) ]  )  
  w = w [ names(w)!='' ]
  w = w / sum( w ) 
  
  nu <- ndesc[u] 
  tu =  descendantSids[[u]] 
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


#' external wrapper for inner_compare_age_groups
#'
#' @export 
compare_age_groups <- function( targetnodes=NULL , scanner_env=readRDS("scanner-env-2021-03-03.rds"), 
                                path_to_data = '/cephfs/covid/bham/results/phylogenetics/latest/' , 
                                include_pillar1=F, min_date = NULL, max_date = NULL,
                                threshold_growth = .5,fast_return=F) 
{
  
  {## all this could be (pre)computed outside this function
  
  e1 = as.environment( scanner_env )
  attach( e1 )
  
  n <- Ntip( tre ) 
  nnode = Nnode( tre )
  
  ## needed for get_comparator_sample
  ##!! setting to environment
  ndesc <<- sapply( 1:(n+nnode), function(u) length( descendantSids[[u]] ) )
  descsts <<- lapply( 1:(n+nnode), function(u) sts[ na.omit( descendantSids[[u]] )  ]  )
  
  ## regen for get_comparator_sample
  amd <- read.csv( list.files(  paste0( path_to_data , '/alignments/') , patt = 'cog_[0-9\\-]+_metadata.csv', full.names=TRUE) 
                   , stringsAs=FALSE )
  amd$sample_time = decimal_date ( as.Date( amd$sample_date ))
  
  # exclude global 
  amd <- amd[ amd$country == 'UK', ]
  # exclude p1 
  if ( !include_pillar1 )
  {
    lhls <- paste0( '.*/', c( 'MILK', 'ALDP', 'QEUH', 'CAMC'), '.*')
    lhpatt = paste( lhls, collapse = '|' )
    amd <- amd[ grepl(amd$sequence_name, patt = lhpatt ) , ] 
  }
  
  max_time <- Inf 
  if ( !is.null( max_date )){
    max_time <- decimal_date( max_date )
  } else{
    max_date <- Sys.Date()
  }
  
  min_time <- -Inf 
  if (!is.null( min_date ))
    min_time <- decimal_date( min_date )
  
  # sample time 
  ## needed for get_comparator_sample
  ##!! setting to environment
  amd$sts <- decimal_date ( as.Date( amd$sample_date ) )
  amd <- amd [ (amd$sts >= min_time) & (amd$sts <= max_time) , ] 
  amd <<- amd[ !is.na( amd$sequence_name ) , ]
  }##
  
  ## select nodes
  if(is.null(targetnodes)) targetnodes <- inner_condense_clusters(Y, threshold_growth)
  
  ## inner function call
  ret <- list()
  for(u in 1:length(targetnodes))
    ret[[u]] <- inner_compare_age_groups(targetnodes[u], fast_return)
  ##
  detach( e1 )
  return(ret)
}

#' internal function calls only, assuming amd is already in the environment
inner_compare_age_groups <- function(u=406318, fast_return=F){
  
  tu =  descendantSids[[u]] 
  stu = sts[ na.omit( descendantSids[[u]] )  ]
  minstu = min(na.omit(stu ))
  maxstu = max(na.omit(stu) )
  
  comparator_sample <- inner_get_comparator_sample(u, nX = 1000)
  
  ## covert cog ages to phe age bands
  ##!! make up some ages
  age_options <- c("0-4","'5-9","'10-19","20-29","30-39","40-49","50-59","60-69","70-79","80+")
  amd$patient_age <- sample(age_options,size=nrow(amd),replace=T)
  amd1 = amd[ (amd$sample_time >= minstu) & (amd$sample_time <= maxstu), c('sequence_name', 'patient_age') ]
  amd1 <- subset(amd1,!is.na(patient_age)&patient_age!='')
  amd1$age <- match(amd1$patient_age,age_options)
  # agemap <- rep(1:10,times=c(5,5,10,10,10,10,10,10,10,45))
  # amd1$age <- agemap[1+amd1$source_age]
  age_samples <- amd1$age[amd1$sequence_name%in%tu]
  n_samples <- length(age_samples)
  comparator_age_samples <- amd1$age[amd1$sequence_name%in%comparator_sample]
  n_comp <- length(comparator_age_samples)
  
  ## test statistic
  all_samples <- c(age_samples,comparator_age_samples)
  n_all <- n_samples + n_comp
  cumulative1 <- sapply(1:10,function(x)sum(age_samples<=x)/n_samples)
  cumulative2 <- sapply(1:10,function(x)sum(all_samples<=x)/n_all)#
  s1stat <- max(abs(cumulative1 - cumulative2))
  
  ## fast return: Birnbaum et al 1951
  # https://doi.org/10.1214/aoms/1177729550
  if(fast_return) {
    x <- s1stat
    n <- n_samples
    jj <- seq.int(from = 0, to = floor(n * (1 - x)))
    osss <- 1 - x * sum(exp(lchoose(n, jj) + (n - jj) * log(1 - x - jj/n) + (jj - 1) * log(x + jj/n)))
    pval <- 1 - osss
    return(pval)
  }
  
  # probability distributions / resampling
  # smpl <- 480; 
  # pad2 <- factor(c(rep(0,sum(all_samples==1)),all_samples,rep(11,sum(all_samples==10))))
  # x2 <- smicd::kdeAlgo(pad2,classes=as.numeric(c(levels(pad2),12)),evalpoints = 100)
  # probs2 <- approx(x=x2$gridx,y=x2$resultDensity[,smpl],xout=2:11)$y
  # cumulative2 <- cumsum(probs2/sum(probs2))
  nboot <- 1000
  nullstat <- c()
  for(i in 1:nboot){
    s1froms2 <- sample(x=all_samples,size=n_samples,replace = T)
    #s1froms2 <- sample(x=1:10,size=n_samples,prob=probs2,replace = T)
    cumulative1 <- sapply(1:10,function(x)sum(s1froms2<=x)/n_samples)
    nullstat[i] <- max(abs(cumulative1 - cumulative2))
  }
  pval <- sum(nullstat>s1stat)/nboot
  return(pval)
  
}

#'
#'
#' @export 
get_mlesky <- function( u=406318 , scanner_env=readRDS("scanner-env-2021-03-03.rds")) 
{
  require(mlesky)
  
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
  
  msg
}
