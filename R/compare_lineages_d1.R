#  _____ ___    ____  _   _ _   _    ___  _   _    ____ _     ___ __  __ ____  
# |_   _/ _ \  |  _ \| | | | \ | |  / _ \| \ | |  / ___| |   |_ _|  \/  | __ ) 
#   | || | | | | |_) | | | |  \| | | | | |  \| | | |   | |    | || |\/| |  _ \ 
#   | || |_| | |  _ <| |_| | |\  | | |_| | |\  | | |___| |___ | || |  | | |_) |
#   |_| \___/  |_| \_\\___/|_| \_|  \___/|_| \_|  \____|_____|___|_|  |_|____/ 



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
  td = datetree( mltr , civmd = civmd, meanrate = meanrate ) 
  res = diff( range( epiweek( date_decimal( td$sts )  )  )  ) + 1
  res <- res * 2
  sg = mlskygrid( td, tau = NULL, tau_lower=.001, tau_upper = 10 , sampleTimes = td$sts , res = res, ncpu = 6, NeStartTimeBeforePresent = 0.25)
  res = with( sg, approx( time, ne, rule = 1, xout = taxis )$y )
  print( res ) 
  res 
}



#' Comparing two lineages in terms of Ne and R using MLESKY. Uses treedater to bootstrap.
#' #'
#' @param trefn main lineage to plot
#' @param matchedfn comparison lineage to plot 
#' @param meanrate mean clock rate for treedater
#' @param taxis time bounds of the analysis
dater_mlesky_plot <- function(trefn, matchedfn, ncpu = 5, meanrate = .001,   taxis = decimal_date( seq( as.Date( '2020-10-15') , as.Date('2021-01-03'), by = 1) )) {
  
  lineage_main = lapply(strsplit(tail(strsplit(trefn, '/')[[1]], 1), '_'), function(x) { paste(toString(x[[1]]), toString(x[[2]]), sep = '_') })[[1]]
  lineage_matched = lapply(strsplit(tail(strsplit(matchedfn, '/')[[1]], 1), '_'), function(x) { paste(toString(x[[1]]), toString(x[[2]]), sep = '_') })[[1]]
  
  
  ofn <- glue::glue( 'd1-tN-{lineage_main}_vs_{lineage_matched}_meanrate{meanrate}{ifelse( grepl(trefn,patt="dedup") , "-dedup", "")}.rds'  )
  dedup <- glue::glue( '{ifelse( grepl(trefn,patt="dedup") , "-dedup", "")}')
  
  tres <- read.tree( trefn ) 
  mtres <- read.tree( matchedfn )
  
  # metadata 
  civetfn =  list.files(  '/cephfs/covid/bham/climb-covid19-volze/phylolatest/civet/' , patt = 'cog_global_[0-9\\-]+_metadata.csv', full.names=TRUE) #'../phylolatest/civet/cog_global_2020-12-01_metadata.csv'
  civmd = read.csv( civetfn , stringsAs=FALSE , header=TRUE )
  civmd$central_sample_id <-  sapply( strsplit( civmd$sequence_name , split='/' ) , '[', 2 ) # for linkage 
  civmd$sample_date <- as.Date( civmd$sample_date )
  civmd$sample_time <- decimal_date( civmd$sample_date ) 
  
  # checking all samples have metadata attached... removing tips that aren't able to be matched
  tres = lapply(tres, function(tr) ape::drop.tip(tr, tr$tip.label[!tr$tip.label %in% civmd$central_sample_id]))
  mtres = lapply(mtres, function(tr) ape::drop.tip(tr, tr$tip.label[!tr$tip.label %in% civmd$central_sample_id]))
  
  
  nes <- tryCatch( {  parallel::mclapply( tres, bootrep, civmd = civmd, meanrate = meanrate, taxis = taxis, mc.cores =  ncpu ) },
                   error = function(e){
                     print(nes)
                     stop(' \n bootrep didnt work... \n')
                   }) 
  N <- do.call( cbind, nes ) 
  
  m_nes <- tryCatch( {  parallel::mclapply( mtres, bootrep, civmd = civmd, meanrate = meanrate,  taxis = taxis, mc.cores =  ncpu ) },
                     error = function(e){
                       print(mtres)
                       stop(' \n bootrep didnt work... \n')
                     })
  m_N <- do.call( cbind, m_nes ) 
  
  saveRDS( list( time = taxis, ne = N, mane = m_N ) , file=ofn )
  
  return(list(pldf = pldf, 
              gpldf = gpldf,
              Rrpldf = Rrpldf,
              q_gr = q_gr,
              q_magr = q_magr,
              qRr = qRr))
  
  
}



# #############
#  ____   ___ ____  _        ___  _       _ ____  
# |___ \ / _ \___ \/ |      / _ \/ |     / | ___| 
#   __) | | | |__) | |_____| | | | |_____| |___ \ 
#  / __/| |_| / __/| |_____| |_| | |_____| |___) |
# |_____|\___/_____|_|      \___/|_|     |_|____/ 


# These calls will save an RDS like "d1-tN-sampler1_B.1.1.7_vs_matchSample_B.1.177_meanrate0.001.rds"

# meanrate = .001
B.1.1.7_matched_notB.1.1.7_meanrate_0.001 = dater_mlesky_plot(trefn = '/cephfs/covid/bham/climb-covid19-boydo/coguk_subsamples/B.1.1.7/sampler1_B.1.1.7_2021-01-15_n=1000.nwk',
                                                              matchedfn = '/cephfs/covid/bham/climb-covid19-boydo/coguk_subsamples/notB.1.1.7/matchSample_notB.1.1.7_2021-01-15_n=1000.nwk', ncpu = 5, meanrate = .001)




B.1.1.7_matched_notB.1.1.7_deduped_meanrate_0.001 = dater_mlesky_plot(trefn = '/cephfs/covid/bham/climb-covid19-boydo/coguk_subsamples/B.1.1.7/sampler1_B.1.1.7_2021-01-15_n=1000-deduped.nwk',
                                                                      matchedfn = '/cephfs/covid/bham/climb-covid19-boydo/coguk_subsamples/notB.1.1.7/matchSample_notB.1.1.7_2021-01-15_n=1000-deduped.nwk', ncpu = 5, meanrate = .001)


#!
B.1.1.7_B.1.177_meanrate_0.001 = dater_mlesky_plot(trefn = '/cephfs/covid/bham/climb-covid19-boydo/coguk_subsamples/B.1.1.7/sampler1_B.1.1.7_2021-01-15_n=1000.nwk',
                                                   matchedfn = '/cephfs/covid/bham/climb-covid19-boydo/coguk_subsamples/B.1.177/matchSample_B.1.177_2021-01-15_n=1000.nwk', ncpu = 5, meanrate = .001)


B.1.1.7_B.1.177_deduped_meanrate_0.001 = dater_mlesky_plot(trefn = '/cephfs/covid/bham/climb-covid19-boydo/coguk_subsamples/B.1.1.7/sampler1_B.1.1.7_2021-01-15_n=1000-deduped.nwk',
                                                           matchedfn = '/cephfs/covid/bham/climb-covid19-boydo/coguk_subsamples/B.1.177/matchSample_B.1.177_2021-01-15_n=1000-deduped.nwk', ncpu = 5, meanrate = .001)



# meanrate = .0005

B.1.1.7_matched_notB.1.1.7_meanrate_0.0005 = dater_mlesky_plot(trefn = '/cephfs/covid/bham/climb-covid19-boydo/coguk_subsamples/B.1.1.7/sampler1_B.1.1.7_2021-01-15_n=1000.nwk',
                                                               matchedfn = '/cephfs/covid/bham/climb-covid19-boydo/coguk_subsamples/notB.1.1.7/matchSample_notB.1.1.7_2021-01-15_n=1000.nwk', ncpu = 5,
                                                               meanrate = .0005)




B.1.1.7_matched_notB.1.1.7_deduped_0.0005 = dater_mlesky_plot(trefn = '/cephfs/covid/bham/climb-covid19-boydo/coguk_subsamples/B.1.1.7/sampler1_B.1.1.7_2021-01-15_n=1000-deduped.nwk',
                                                              matchedfn = '/cephfs/covid/bham/climb-covid19-boydo/coguk_subsamples/notB.1.1.7/matchSample_notB.1.1.7_2021-01-15_n=1000-deduped.nwk', ncpu = 5,
                                                              meanrate = .0005)



#!
B.1.1.7_B.1.177_0.0005 = dater_mlesky_plot(trefn = '/cephfs/covid/bham/climb-covid19-boydo/coguk_subsamples/B.1.1.7/sampler1_B.1.1.7_2021-01-15_n=1000.nwk',
                                           matchedfn = '/cephfs/covid/bham/climb-covid19-boydo/coguk_subsamples/B.1.177/matchSample_B.1.177_2021-01-15_n=1000.nwk', ncpu = 5,
                                           meanrate = .0005)


B.1.1.7_B.1.177_deduped_0.0005 = dater_mlesky_plot(trefn = '/cephfs/covid/bham/climb-covid19-boydo/coguk_subsamples/B.1.1.7/sampler1_B.1.1.7_2021-01-15_n=1000-deduped.nwk',
                                                   matchedfn = '/cephfs/covid/bham/climb-covid19-boydo/coguk_subsamples/B.1.177/matchSample_B.1.177_2021-01-15_n=1000-deduped.nwk', ncpu = 5,
                                                   meanrate = .0005)





