
datetree <- function(mltr, civmd, meanrate)
{
  sts <- setNames( civmd$sample_time[  match( mltr$tip.label , civmd$central_sample_id ) ], mltr$tip.label )
  tr <- di2multi(mltr, tol = 1e-05)
  tr = unroot(multi2di(tr))
  tr$edge.length <- pmax(1/29000/5, tr$edge.length)
  dater(unroot(tr), sts[tr$tip.label], s = 29000, omega0 = meanrate, numStartConditions = 0, meanRateLimits = c(meanrate, meanrate + 1e-6), ncpu = 6)
}



bootrep <- function(mltr, civmd = civmd, meanrate = meanrate, taxis = taxis, ...)
{
  td = datetree( mltr , civmd = civmd, meanrate = meanrate, ... ) 
  res = diff( range( epiweek( date_decimal( td$sts )  )  )  ) + 1
  res <- res * 2
  sg = mlskygrid( td, tau = NULL, tau_lower=.001, tau_upper = 10 , sampleTimes = td$sts , res = 12, ncpu = 6)
  res = with( sg, approx( time, ne, rule = 1, xout = taxis )$y )
  print( res ) 
  res 
}



#' Writes bash scripts to be launched on the CLIMB server via slurm scheduler. Scripts will recomputing ML tree using IQtree from an alignment (expressed in the form of a newick tree)
#'
#' @param tres tree whose tips are to form the ML tree
recomputing_ML_trees = function(tres = '/cephfs/covid/bham/climb-covid19-volze/subsampling/sampler1_B.1.1.7_2021-01-04_n=1000.nwk',
                                root_dir = '..') {
  
  
  
  tmp = tail(strsplit(tres, '/')[[1]], 1)
  

  lineage = lapply(strsplit(tmp, '_'), function(x) { paste(toString(x[[1]]), toString(x[[2]]), sep = '_') })[[1]]
  
  
  
  # READING IN GLOBAL ALIGNMENT
  cat(' \n READING IN GLOBAL ALIGNMENT... \n')
  algn_fn <- 	list.files(  paste0(root_dir, '/phylolatest/civet/' ), patt = 'cog_global_[0-9\\-]+_alignment.fasta', full.names=TRUE) #'../phylolatest/civet/cog_global_2020-12-01_metadata.csv'
  t0_readalgn = Sys.time()
  algn <- read.dna(algn_fn, format = 'fasta' ) 
  t1_readalgn = Sys.time()
  t1_readalgn-t0_readalgn
  
  
  
  
  # READING IN ERIKS TREE WHOSE SAMPLES I MATCH
  X = read.tree(tres)
  
  ncpu=24
  
  
  Xs <- lapply(X, function(x) {
    return(list(central_sample_id = x$tip.label))
  } )
  
  for(i in 1:length(Xs)) {
    Xs[[i]]$replicate = i
    Xs[[i]]$sample_size = 1000
  }

  # 
  cat(' \n WRITING FASTAS... \n')
  
  t0 = Sys.time()
  parallel::mclapply(Xs, function(i ) {
    write.dna(algn[lapply(strsplit(rownames(algn), '/'), '[[', 2) %in% as.character(i$central_sample_id), ], file = 
                paste0(lineage, '_samplesize', i$sample_size[1], 'rep', i$replicate[1], '.fasta'), format = 'fasta')}, 
    mc.cores = ncpu)
  
  
  bash_template = 	readLines( '~/submit.sh' ) 
  
  cat(' \n WRITING bash scripts \n')
  
  # writing bash scripts
  write_bash_for_iqtree <- function(x,  bash_template, lineage) {
    
    bashofn = gsub( bash_template, pattern='test_apply', replacement = paste0('job_', lineage, '_samplesize', x$sample_size, 'rep', x$replicate) )
    bashofn = gsub( bashofn, pattern='<algn>', replacement = paste0(getwd(), '/' , lineage, '_samplesize', x$sample_size, 'rep', x$replicate, '.fasta'))
    writeLines(bashofn, con = paste0('job_', lineage, '_samplesize', x$sample_size, 'rep', x$replicate, '.sh'))
    
  }

  lapply(Xs, write_bash_for_iqtree,  bash_template = bash_template, lineage = lineage)
  
  
  return(list(Xs = Xs, lineage = lineage))
  
}


#' recomputing ML tree using IQtree from an alignment
#' #'
#' @param Xs list of jobs returned from recomputing_ML_trees
running_IQ_tree_slurm <- function(Xs, lineage) {
  
  # launching IQTREE jobs
  lapply(Xs, function(x) {
    system(paste0('sbatch job_', lineage, '_samplesize', x$sample_size, 'rep', x$replicate, '.sh'))
  })
  
}



#' Comparing two lineages in terms of Ne and R using MLESKY. Uses treedater to bootstrap.
#' #'
#' @param trefn main lineage to plot
#' @param matchedfn comparison lineage to plot (can be matched or unmatched)
dater_mlesky_plot <- function(trefn, matchedfn, ncpu = 5, meanrate = .001) {
  
  
  
  
  lineage_main = lapply(strsplit(tail(strsplit(trefn, '/')[[1]], 1), '_'), function(x) { paste(toString(x[[1]]), toString(x[[2]]), sep = '_') })[[1]]
  lineage_matched = lapply(strsplit(tail(strsplit(matchedfn, '/')[[1]], 1), '_'), function(x) { paste(toString(x[[1]]), toString(x[[2]]), sep = '_') })[[1]]
  
 
  ofn <- glue::glue( 'd1-tN-{lineage_main}_vs_{lineage_matched}_meanrate{meanrate}{ifelse( grepl(trefn,patt="dedup") , "-dedup", "")}.rds'  )
  
  dedup <- glue::glue( '{ifelse( grepl(trefn,patt="dedup") , "-dedup", "")}')
  
  tres <- read.tree( trefn ) 
  mtres = read.tree( matchedfn )
  

  # put results on common time axis 
  taxis <- decimal_date( seq( as.Date( '2020-10-15') , as.Date('2020-12-18')  , by = 1) )
  
  # metadata 
  civetfn =  list.files(  '/cephfs/covid/bham/climb-covid19-volze/phylolatest/civet/' , patt = 'cog_global_[0-9\\-]+_metadata.csv', full.names=TRUE) #'../phylolatest/civet/cog_global_2020-12-01_metadata.csv'
  civmd = read.csv( civetfn , stringsAs=FALSE , header=TRUE )
  
  
  civmd$central_sample_id <-  sapply( strsplit( civmd$sequence_name , split='/' ) , '[', 2 ) # for linkage 
  civmd$sample_date <- as.Date( civmd$sample_date )
  civmd$sample_time <- decimal_date( civmd$sample_date ) 
  
  # checking all samples have metadata attached... removing tips that aren't able to be matched
  tres = lapply(tres, function(tr) ape::drop.tip(tr, tr$tip.label[!tr$tip.label %in% civmd$central_sample_id]))
  mtres = lapply(mtres, function(tr) ape::drop.tip(tr, tr$tip.label[!tr$tip.label %in% civmd$central_sample_id]))
  

  nes <- tryCatch( {  parallel::mclapply( tres, bootrep, civmd = civmd, meanrate = meanrate, taxis = taxis, mc.cores =  3 ) },
                   error = function(e){
                     print(nes)
                     stop(' \n bootrep didnt work... \n')
                   }
                   
  ) 
  
  
  N <- do.call( cbind, nes ) 
  
  m_nes <- tryCatch( {  parallel::mclapply( mtres, bootrep, civmd = civmd, meanrate = meanrate,  taxis = taxis, mc.cores =  3 ) },
                     error = function(e){
                       print(mtres)
                       stop(' \n bootrep didnt work... \n')
                     }
                     
  )
  
  
  m_N <- do.call( cbind, m_nes ) 
  
  saveRDS( list( time = taxis, ne = N, mane = m_N )
           , file=ofn )
  

  
  print(
    quantile( apply( N, MAR=2, FUN = function(x) median( exp(diff(log(na.omit(x)))) )^6.5  ), c(.025, .5, .975)) 
  )

  
  #plot results 
# 
#    _   _  ___ _____ _____ 
#   | \ | |/ _ \_   _| ____|
#   |  \| | | | || | |  _|  
#   | |\  | |_| || | | |___ 
#   |_| \_|\___/ |_| |_____|
#     
#   better plots are found in N501_analyses/analyse_results_two_lineages_after_bootstrapping_and_mlesky.R
  
  
  tN = readRDS( ofn )
  # tN = readRDS( 'regenerated_ML_trees_d1-tN-meanrate0.001.rds' )
  attach( tN )
  q_ne = t(apply( ne, 1, function(x) quantile( na.omit(x), c(.5, .025, .975 )) ))
  q_mane  = t(apply( mane, 1, function(x) quantile( na.omit(x), c(.5, .025, .975 )) ))
  
  gr =  apply( ne, 2, function(x) c(exp( diff(log(x)) )^6.5, NA)  ) 
  magr =  apply( mane, 2, function(x) c(exp( diff(log(x)) )^6.5, NA)  ) 
  
  q_gr = t( apply( gr, 1, function(x) quantile( na.omit(x), c(.5, .025, .975 )) ) )
  q_magr  = t( apply( magr, 1, function(x) quantile( na.omit(x), c(.5, .025, .975 )) ) )
  
  colnames( q_ne ) = colnames( q_mane ) = c( 'y', 'ylb', 'yub' )
  pldf0 = as.data.frame( q_ne ) ; pldf0$lineage = lineage_main; pldf0$time = time 
  pldf1 = as.data.frame( q_mane ); pldf1$lineage = lineage_matched; pldf1$time = time 
  pldf = rbind( pldf0, pldf1 )
  
  p0 = ggplot( aes(x = as.Date( date_decimal( time)), y = y, colour = lineage, fill = lineage , ymin = ylb, ymax = yub ) , data = pldf ) + geom_path() + geom_ribbon( alpha = .25 ) + xlab('') + ylab('Effective population size' ) + theme_minimal() + theme(legend.position='none')
  
  colnames( q_gr ) = colnames( q_magr ) = c( 'y', 'ylb', 'yub' )
  gpldf0 = as.data.frame( q_gr ) ; gpldf0$lineage = lineage_main; gpldf0$time = time 
  gpldf1 = as.data.frame( q_magr ); gpldf1$lineage = lineage_matched; gpldf1$time = time 
  gpldf = rbind( gpldf0, gpldf1 )
  p1 =  ggplot( aes(x = as.Date( date_decimal( time)), y = y, colour = lineage, fill = lineage , ymin = ylb, ymax = yub ) , data = gpldf ) + geom_path() + geom_ribbon( alpha = .25 ) + xlab('') + ylab('Reproduction number' ) + theme_minimal() + scale_y_log10() + geom_hline( aes( yintercept=1), lty = 2 )  + theme(legend.position='top')
  
  Rratio = gr / magr 
  qRr =  t( apply( Rratio, 1, function(x) quantile( na.omit(x), c(.5, .025, .975 )) ) )
  colnames( qRr ) =  c( 'y', 'ylb', 'yub' )
  Rrpldf = as.data.frame( qRr ); Rrpldf$time = time 
  p2 = ggplot( aes(x = as.Date( date_decimal( time)), y = y, ymin = ylb, ymax = yub ) , data = Rrpldf ) + geom_path() + geom_ribbon( alpha = .25 ) + xlab('') + ylab('Ratio of reproduction numbers' ) + theme_minimal() + scale_y_log10() + geom_hline( aes( yintercept=1), lty = 2 ) 
  
  library( cowplot )
  P0 = plot_grid( plotlist = list( p0, p1, p2 ), nrow = 1 )
  # ggsave( plot = P0, file = 'd1.png', width = 8, height = 3 )
  ggsave( plot = P0, file = paste0('d1_', lineage_main, '_', lineage_matched, dedup, '.pdf'), width = 8, height = 3 )
  # ggsave( plot = P0, file = 'd1_regen.pdf', width = 8, height = 3 )
  
  
  return(list(pldf = pldf, 
              gpldf = gpldf,
              Rrpldf = Rrpldf,
              q_gr = q_gr,
              q_magr = q_magr,
              qRr = qRr))
  
  
}


