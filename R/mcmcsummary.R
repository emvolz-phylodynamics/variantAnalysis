# code adapted from bayesiantools package 
# https://github.com/florianhartig/BayesianTools/R/classMcmcSamplerList.R

summary.mcmc <- function(object, ...){
  #codaChain = getSample(sampler, parametersOnly = parametersOnly, coda = T, ...)
  #summary(codaChain)
  #rejectionRate(sampler$codaChain)
  #effectiveSize(sampler$codaChain)
  #DIC(sampler)
  #max()
  #}
  
  sampler <- object
  
  try(DInf <- DIC(sampler), silent = TRUE)
  MAPvals <- round(MAP(sampler)$parametersMAP,3)
  psf <- FALSE
  
  mcmcsampler <- sampler$settings$sampler
  runtime <- sampler$settings$runtime[3]
  correlations <- round(cor(getSample(sampler)),3)
  
  chain <- getSample(sampler, parametersOnly = T, coda = T, ...)
  # chain <- getSample(sampler, parametersOnly = T, coda = T)
  if("mcmc.list" %in% class(chain)){
    psf <- TRUE
    nrChain <- length(chain)
    nrIter <- nrow(chain[[1]])
    conv <- ifelse(chain$setup$numPars > 1, round(coda::gelman.diag(chain)$mpsrf,3), round(coda::gelman.diag(chain)$mpsrf,3)$psrf[1])
    npar <- sampler$setup$numPars
    lowerq <- upperq <- numeric(npar)
    medi <- numeric(npar)
    parnames <- colnames(chain[[1]])
    
    # Shorthen parameter names
    for (i in 1:npar) {
      if (nchar(parnames[i]) > 8)
        parnames[i] <- paste(substring(parnames[i], 1, 6), "...", sep = "")
    }
    
    for (i in 1:npar) {
      tmp <- unlist(chain[, i])
      tmp <- quantile(tmp, probs = c(0.025, 0.5, 0.975))
      lowerq[i] <- round(tmp[1], 3)
      medi[i] <- round(tmp[2], 3)
      upperq[i] <- round(tmp[3], 3)
    }
    
  } else{
    nrChain <- 1
    nrIter <- nrow(chain)
    npar <- sampler$setup$numPars
    conv <- "Only one chain; convergence cannot be determined!"
    medi <- numeric(npar)
    lowerq <- upperq <- numeric(npar)
    parnames <- colnames(chain)
    for (i in 1:npar) {
      tmp <- quantile(chain[, i], probs = c(0.025, 0.5, 0.975))
      lowerq[i] <- round(tmp[1], 3)
      medi[i] <- round(tmp[2], 3)
      upperq[i] <- round(tmp[3], 3)
    }
    
  }
  
  parOutDF <- cbind(MAPvals, lowerq, medi, upperq)
  colnames(parOutDF) <- c("MAP", "2.5%", "median", "97.5%")
  if (psf == TRUE) {
    psf <- round(gelmanDiagnostics(sampler)$psrf[,1], 3)
    parOutDF <- cbind(psf, parOutDF)
  }
  row.names(parOutDF) <- parnames
  
  cat(rep("#", 25), "\n")
  cat("## MCMC chain summary ##","\n")
  cat(rep("#", 25), "\n", "\n")
  cat("# MCMC sampler: ",mcmcsampler, "\n")
  cat("# Nr. Chains: ", nrChain, "\n")
  cat("# Iterations per chain: ", nrIter, "\n")
  cat("# Rejection rate: ", ifelse(sampler$setup$numPars == 1, round(mean(sapply(chain, coda::rejectionRate)),3), round(mean(coda::rejectionRate(chain)),3) ) , "\n")
  cat("# Effective sample size: ", ifelse(sampler$setup$numPars == 1, round(coda::effectiveSize(chain),0), round(mean(coda::effectiveSize(chain)),0) ) , "\n")
  cat("# Runtime: ", runtime, " sec.","\n", "\n")
  cat("# Parameters\n")
  print(parOutDF)
  cat("\n")
  
  try(cat("## DIC: ", round(DInf$DIC,3), "\n"), silent = TRUE)
  cat("## Convergence" ,"\n", "Gelman Rubin multivariate psrf: ", conv, "\n","\n")
  invisible( object )
}

# for list 
summary.mcmclist <- function(object, ...){
  #codaChain = getSample(sampler, parametersOnly = parametersOnly, coda = T, ...)
  #summary(codaChain)
  #rejectionRate(sampler$codaChain)
  #effectiveSize(sampler$codaChain)
  #DIC(sampler)
  #max()

  sampler <- object

  DInf <- DIC(sampler)
  MAPvals <- round(MAP(sampler)$parametersMAP,3)

  gelDiag <- gelmanDiagnostics(sampler)
  psf <- round(gelDiag$psrf[,1], 3)
  
  mcmcsampler <- sampler[[1]]$settings$sampler
  
  runtime <- 0
  for(i in 1:length(sampler)) runtime <- runtime+sampler[[i]]$settings$runtime[3]

  correlations <- round(cor(getSample(sampler)),3)

  
  sampler <- getSample(sampler, parametersOnly = T, coda = T)
  if("mcmc.list" %in% class(sampler)){
    nrChain <- length(sampler)
    nrIter <- nrow(sampler[[1]])
    conv <- round(gelDiag$mpsrf,3)
    npar <- ncol(sampler[[1]])
    lowerq <- upperq <- numeric(npar)
    medi <- numeric(npar)
    parnames <- colnames(sampler[[1]])
    for(i in 1:npar){
      tmp <- unlist(sampler[,i])
      tmp <- quantile(tmp, probs = c(0.025, 0.5, 0.975))
      lowerq[i] <- round(tmp[1],3)
      medi[i] <- round(tmp[2],3)
      upperq[i] <- round(tmp[3],3)
    }

  }else{
    nrChain <- 1
    nrIter <- nrow(sampler)
    npar <- ncol(sampler)
    conv <- "Only one chain; convergence cannot be determined!"
    medi <- numeric(npar)
    lowerq <- upperq <- numeric(npar)
    parnames <- colnames(sampler)
    for(i in 1:npar){
      tmp <- quantile(sampler[,i], probs = c(0.025, 0.5, 0.975))
      lowerq[i] <- round(tmp[1],3)
      medi[i] <- round(tmp[2],3)
      upperq[i] <- round(tmp[3],3)
    }

  }
  
  # output for parameter metrics
  parOutDF <- cbind(psf, MAPvals, lowerq, medi, upperq)
  colnames(parOutDF) <- c("psf", "MAP", "2.5%", "median", "97.5%")
  row.names(parOutDF) <- parnames

  
  cat(rep("#", 25), "\n")
  cat("## MCMC chain summary ##","\n")
  cat(rep("#", 25), "\n", "\n")
  cat("# MCMC sampler: ",mcmcsampler, "\n")
  cat("# Nr. Chains: ", nrChain, "\n")
  cat("# Iterations per chain: ", nrIter, "\n")
  cat("# Rejection rate: ", round(mean(coda::rejectionRate(sampler)),3), "\n")
  cat("# Effective sample size: ", round(mean(coda::effectiveSize(sampler)),0), "\n")
  cat("# Runtime: ", runtime, " sec.","\n", "\n")
  cat("# Parameters\n")
  print(parOutDF)
  cat("\n")
  cat("## DIC: ", round(DInf$DIC,3), "\n")
  cat("## Convergence" ,"\n", "Gelman Rubin multivariate psrf: ", conv, "\n","\n")
  
invisible( object )   
}


#' @export 
hier_bayes_exponentialGrowth_frequency.summary2 <- function(f, d= NULL,  MU = 73, thin = 1, start = 1000)
{
	library( BayesianTools )
	library( ggplot2 )
	library( ggforce )
	
	if ( is.null( d ) )
		d = f$d 
	
	print( (fsum = summary.mcmc( f )  ) )
	pnames =  colnames( fsum$X )
	M <- (length( pnames )-1)/2


	

	print( 'cv ' )
	print( quantile( getSample( f$chain[[1]] , start = start, which=pnames=='cv' ) , c( .5, .025, .975 ) ) ) 

	pnames2 <- c( pnames, c( 'LP', 'L', 'Pr' ) )
	lp <- sapply( f$chain, function(chain){
		colnames(chain) <- pnames2 
		getSample( chain[, 'LP'], thin = thin, start = start)
	}) 
	X11(); matplot( lp, type = 'l' , main = 'log Posterior') 
	selcoef2 <-   sapply( f$chain, function(chain){
		colnames(chain) <- pnames2 
		rs <- getSample( chain, thin = thin, start = start )[, (M+1):(M+M)]
		drs <- rowMeans( rs[ , d$genotype=='wt' ] )
		grs <- rowMeans( rs[ , d$genotype=='mutant' ] )
		.selcoef <- ( grs - drs ) / MU 
		.selcoef
	}) 
	X11(); matplot( selcoef2, type = 'l' , main = 'Selection coefficient' ) 
	cat('Selection coefficient \n') 
	print( quantile( selcoef2[,1], c(.5, .025, .975 )) )
	
	#map 
	x = MAP( f ) 
	theta <- x[[1]] 
		logns <- theta[1:M] 
	rs <- theta[ (M+1):(M+M) ]
	drs <- rs[ d$genotype=='wt' ]
	grs <- rs[ d$genotype=='mutant' ]	
	dmu = mean( drs ) 
	gmu = mean( grs ) 
	selcoef <- (gmu - dmu) / MU
	sigma <- theta['cv'] * MU
	r.squared_map <- 1 - (sum( (drs - dmu)^2 )  + sum( (grs - gmu)^2) )/sum((rs - mean(rs))^2)

	#figs
	chain = f$chain[[1]]
	colnames(chain) <- pnames2 
	
	densd <- data.frame( wt = rnorm( 5e4, dmu, sigma ), mutant=rnorm( 5e4,gmu,sigma) )
	p = ggplot(data = densd) + geom_density(aes(x=wt, y = ..density..) , fill='#7d8c85', alpha = .5) + 
		geom_density(aes(x=mutant, y = -..density..) , fill='#c8874e', alpha =.5)
	
	rdf <- data.frame( r = rs, genotype = d$genotype
		, position= runif(nrow(d), .001, .0085) * (ifelse(d$genotype=='mutant', -1, 1) )
		, size = unname( d$n )
		, col = ifelse( d$genotype=='wt', '#718e87', '#dc863b' )
		, first_sample = d$first_sample 
		, stringsAsFactors=FALSE
	)
	pp = p + geom_point( aes( r, position, size = sqrt( size), col=genotype), data = rdf ) +  scale_color_manual(values=c( wt='#718e87', mutant='#dc863b')) + theme_classic()+ xlab('') + ylab('')
	pp = pp + coord_flip()  + theme(axis.line.x=element_blank(), axis.text.x=element_blank() , axis.ticks.x=element_blank() )  + xlab( 'Cluster growth rate (1/year)' ) 
	#ggsave( pp, file='a2res0-growth.png' , width=3.75, height = 3.5 )
	#ggsave( pp, file='a2res0-growth.svg' , width=3.75, height = 3.5 )
		
	pp
}
