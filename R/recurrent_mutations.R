
#' Identify recurrent mutations in a tree and quantify offspring at nodes with mixed descendents 
#' 
#' @param tr1 ape phylo, rooted, may have polytomies 
#' @param tipstate character vector with names corresonding to tip labels; gives genotype 
#' @param sts optional vector of sample times named according to tip label; gives approx age of mixed nodes 
#' @return Data frame with id of 'mixed' nodes that have offspring of multiple genotypes and the number of offspring of each genotype 
#' @export 
recurrent_mutations <- function(tr1, tipstate, sts = NULL )
{
	library( ape ) 
	library( sarscov2 )
	library( lubridate )
	
	tipstate = tipstate[ tr1$tip.label ]
	#TODO exclude tips missing state here 
	if ( is.null( sts )) 
		sts <- setNames( rep( Inf, Ntip(tr1)), tr1$tip.label)
	
	if ( length( setdiff( tr1$tip.label , names(tipstate))) > 0 ) 
		stop( 'missing tip states for some samples in tree' )
	if ( length( setdiff( tr1$tip.label , names(sts))) > 0 ) 
		stop( 'missing sample time for some samples in tree' )
	
	dgtrs = lapply(1:(Ntip(tr1)+Nnode(tr1)), function(a){ #slow 
		as.vector( tr1$edge[ tr1$edge[,1]==a ,2] )
	})
	
	hpnodedata <- list() # index by node number, $genotype vector of offspring number 
	poedges <- tr1$edge[ postorder( tr1 ) , ] 
	nodestate <- rep(NA, Ntip(tr1) + Nnode(tr1) ) # genotype or 'mixed' 
	nodestate[ 1:Ntip(tr1) ] <- tipstate 
	noff <- rep( 0, Ntip(tr1) + Nnode(tr1)) 
	noff[ 1:Ntip(tr1) ] <- 1
	minsamptime <- rep( Inf, Ntip(tr1) + Nnode(tr1 )) #earliest sample within the clade 
	minsamptime[ 1:Ntip(tr1)] <- sts [ tr1$tip.label ] 
	poancs = poedges[,1]
	poancs = poancs[!duplicated(poancs)]
	for ( a in poancs ){
		uvs = dgtrs[[a]] 
		nstates <- nodestate[ uvs ] 
		noffs <- noff[ uvs ]
		tt <- table( nstates )
		ttnm <- tt [ names(tt)!='mixed' ] 
		minsamptime[ a ] <- min ( minsamptime[ uvs ] )
		if ( length( tt ) == 1  ){ # all identical , maybe all mixed 
			nodestate[a] <- nstates[1]
			noff[a] <- sum( noffs )
			
		} else if ( sum( nstates=='mixed' ) == 1 ){ # one mixed offspring 
			if ( length( ttnm ) == 1 ){ # remainder not mixed 
				nodestate[a] <- names( ttnm ) [1]
				noff[a] <- sum( noffs )
			} else{ # other offspring mixed 
				nodestate[a] <- 'mixed' 
				noff[a] <- 0 
			}
		} else if ( sum( nstates=='mixed' ) > 1 ){
			nodestate[a] <- 'mixed' 
			noff[a] <- 0 
		} else{ # recurrant mutation 
			nodestate[a] <- 'mixed'	
			noff[a] <- 0 
			print( ttnm ) 
			print( tt )
			hpnodedata[[as.character(a)]] <- setNames( lapply( names( ttnm ), function(genotype) {
				noffs[ nstates == genotype  ] 
			}) ,  names( ttnm ))
			#browser() 
		}
	}
	X = do.call( rbind, lapply( names(hpnodedata), function(nm){
		x = hpnodedata[[nm]] 
		#y = expand.grid( x) 
		xnms = names( x )[1:2] #only 2 genotypes supported
		xx = as.list(  setNames(  sapply( x, sum ) , xnms ) )
		y = data.frame( xx  )
		y$node = nm 
		y$minsamptime <- minsamptime[ as.numeric(nm) ] 
		y
	}) )
	X
}

