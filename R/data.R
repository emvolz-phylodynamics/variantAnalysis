
#' Compute AA genotype at specified postions in alignment
#'
#" @param algn_fn path to most recent COGUK alignment 
#' @param positions Named list of nt coordinates for given genotypes; each set of coordinates should span a codon
#' @param ofn optional path to csv for output 
#' @return data frame with genotypes 
#' @export 
compute_genotypes <- function(  
	algn_fn = 'alignments/cog_2020-11-11_alignment.fasta'
	, positions = list( 
		S614 = 23402:23404
		, S222 = 22226:22228
		, N220 = 28931:28933
		, ORF10 = 29645:29647
		, Y453F = 22919:22921
	) 
	, ofn = NULL 
){
	library( ape ) 
	library( lubridate )
	# aligment: 
	a1 <- read.dna( algn_fn, format = 'fasta' ) 
	
	pos <- do.call( c, lapply( positions, function(x) x[1:3] ) )
	p1 = trans(a1) 
	#apply( as.character( p1) , 2 , table )
	
	
	genotype = g =  as.character( p1  )
	colnames( g ) <- names( positions  )
	
	# save genotypes for future reference 
	gdf = as.data.frame( g ) 
	gdf$central_sample_id <- sapply( strsplit(rownames(p1) , split='/' ) , '[', 2 ) 
	if ( !is.null( ofn ))
		write.csv( gdf, file = 'genotypes.csv', quote=FALSE , row.names=FALSE)
	gdf
}
