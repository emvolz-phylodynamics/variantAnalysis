require(ggplot2)
require(lubridate)
library(cowplot)
library(grid)
library(gridExtra)

unpack_two_Lineages <- function(ofn, Lineage_main, Lineage_matched, dedup, meanrate) {
  tN = readRDS( ofn )
  
  attach( tN )
  q_ne = t(apply( ne, 1, function(x) quantile( na.omit(x), c(.5, .025, .975 )) ))
  q_mane  = t(apply( mane, 1, function(x) quantile( na.omit(x), c(.5, .025, .975 )) ))
  
  gr =  apply( ne, 2, function(x) c(exp( diff(log(x)) )^6.5, NA)  ) 
  magr =  apply( mane, 2, function(x) c(exp( diff(log(x)) )^6.5, NA)  ) 
  
  q_gr = t( apply( gr, 1, function(x) quantile( na.omit(x), c(.5, .025, .975 )) ) )
  q_magr  = t( apply( magr, 1, function(x) quantile( na.omit(x), c(.5, .025, .975 )) ) )
  
  colnames( q_ne ) = colnames( q_mane ) = c( 'y', 'ylb', 'yub' )
  pldf0 = as.data.frame( q_ne ) ; pldf0$Lineage = Lineage_main; pldf0$time = time 
  pldf1 = as.data.frame( q_mane ); pldf1$Lineage = Lineage_matched; pldf1$time = time 
  pldf = rbind( pldf0, pldf1 )
  
  p0 = ggplot( aes(x = as.Date( date_decimal( time)), y = y, colour = Lineage, fill = Lineage , ymin = ylb, ymax = yub ) , data = pldf ) + 
    geom_path(size=1) + geom_ribbon( alpha = .25, col = NA ) + xlab('') + ylab('Effective population size' ) + 
    theme_minimal() + theme(legend.position='none',panel.grid.minor = element_blank())+
    scale_x_date(date_breaks = "1 month", date_labels = '%b')+theme(axis.text=element_text(size=12),
                                                                    axis.title=element_text(size=14)) + scale_y_log10()+ 
    annotation_logticks(colour = 'grey', short = unit(.05, "cm"), mid = unit(.05, "cm"), long = unit(.05, "cm"))
  
  colnames( q_gr ) = colnames( q_magr ) = c( 'y', 'ylb', 'yub' )
  gpldf0 = as.data.frame( q_gr ) ; gpldf0$Lineage = Lineage_main; gpldf0$time = time 
  gpldf1 = as.data.frame( q_magr ); gpldf1$Lineage = Lineage_matched; gpldf1$time = time 
  gpldf = rbind( gpldf0, gpldf1 )
  
  Rratio = gr / magr 
  qRr =  t( apply( Rratio, 1, function(x) quantile( na.omit(x), c(.5, .025, .975 )) ) )
  colnames( qRr ) =  c( 'y', 'ylb', 'yub' )
  Rrpldf = as.data.frame( qRr ); Rrpldf$time = time 
  
  
  r_range = range(c(gpldf$yub, Rrpldf$yub, gpldf$ylb, Rrpldf$ylb), na.rm = T)
  
  p1 =  ggplot( aes(x = as.Date( date_decimal( time)), y = y, colour = Lineage, fill = Lineage , ymin = ylb, ymax = yub ) , data = gpldf ) + 
    geom_path(size=1) + geom_ribbon( alpha = .25 , col = NA) + xlab('') + ylab('Reproduction number' ) + theme_minimal() + 
    scale_y_log10(limits = r_range) + annotation_logticks(colour = 'grey', short = unit(.05, "cm"), mid = unit(.05, "cm"), long = unit(.05, "cm"))   + 
    geom_hline( aes( yintercept=1), lty = 2 )  + theme(legend.position='', panel.grid.minor = element_blank())+
    scale_x_date(date_breaks = "1 month", date_labels = '%b')+
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14))
  
  p2 = ggplot( aes(x = as.Date( date_decimal( time)), y = y, ymin = ylb, ymax = yub ) , data = Rrpldf ) + geom_path(size=1) + 
    
    geom_ribbon( alpha = .25, col = NA ) + xlab('') + ylab('Ratio of reproduction numbers' ) + theme_minimal() + 
    scale_y_log10(limits = r_range) + annotation_logticks(colour = 'grey', short = unit(.05, "cm"), mid = unit(.05, "cm"), long = unit(.05, "cm")) + theme( panel.grid.minor = element_blank())+
    geom_hline( aes( yintercept=1), lty = 2 ) +scale_x_date(date_breaks = "1 month", date_labels = '%b')+
    theme(axis.text=element_text(size=12),  axis.title=element_text(size=14))
  
  legend <- cowplot::get_legend(
    p1 + theme(legend.box.margin = margin(0, 0, 0, 12), legend.position = "top", legend.title = element_text(size=14), legend.text =element_text(size=12) ) 
  )
  
  
  
  P0 = cowplot::plot_grid( plotlist = list( p0, p1, p2 ), nrow = 1, align = "v" )
  
  P0 = cowplot::plot_grid(legend, P0, ncol = 1, rel_heights =  c(0.1, 1))
  
  ggsave( plot = P0, file = paste0('d1_', Lineage_main, '_', Lineage_matched, dedup, '_', "meanrate", meanrate, '.pdf'), width = 12, height = 4.5 )
  
  Time=time
  
  detach(tN)
  return(list(ofn = ofn, pldf = pldf, 
              time = Time,
              gpldf = gpldf,
              Rrpldf = Rrpldf,
              q_gr = q_gr,
              q_magr = q_magr,
              qRr = qRr))
}



# mean rate = 0.001


B.1.1.7_matched_notB.1.1.7_meanrate0.001 = 
  unpack_two_Lineages(ofn = "d1-tN-sampler1_B.1.1.7_vs_matchSample_notB.1.1.7_meanrate0.001.rds",
  Lineage_main = 'B.1.1.7',
  Lineage_matched = 'Control',
  dedup = "",
  meanrate = 0.001)





B.1.1.7_matched_notB.1.1.7_deduped_meanrate0.001 = 
  unpack_two_Lineages(ofn = "d1-tN-sampler1_B.1.1.7_vs_matchSample_notB.1.1.7_meanrate0.001-dedup.rds",
                      Lineage_main = 'B.1.1.7',
                      Lineage_matched = 'Control',
                      dedup = "-dedup",
                      meanrate = 0.001)



B.1.1.7_B.1.177_meanrate0.001 = 
  unpack_two_Lineages(ofn = "d1-tN-sampler1_B.1.1.7_vs_matchSample_B.1.177_meanrate0.001.rds",
                      Lineage_main = 'B.1.1.7',
                      Lineage_matched = 'B.1.177',
                      dedup = "",
                      meanrate = 0.001)



B.1.1.7_B.1.177_deduped_meanrate0.001 = 
  unpack_two_Lineages(ofn = "d1-tN-sampler1_B.1.1.7_vs_matchSample_B.1.177_meanrate0.001-dedup.rds",
                      Lineage_main = 'B.1.1.7',
                      Lineage_matched = 'B.1.177',
                      dedup = "-dedup",
                      meanrate = 0.001)



# mean rate = 0.0005

B.1.1.7_matched_notB.1.1.7_meanrate0.0005 = 
  unpack_two_Lineages(ofn = "d1-tN-sampler1_B.1.1.7_vs_matchSample_notB.1.1.7_meanrate5e-04.rds",
                      Lineage_main = 'B.1.1.7',
                      Lineage_matched = 'Control',
                      dedup = "",
                      meanrate = 0.0005)

B.1.1.7_matched_notB.1.1.7_deduped_meanrate0.0005 = 
  unpack_two_Lineages(ofn = "d1-tN-sampler1_B.1.1.7_vs_matchSample_notB.1.1.7_meanrate5e-04-dedup.rds",
                      Lineage_main = 'B.1.1.7',
                      Lineage_matched = 'Control',
                      dedup = "-dedup",
                      meanrate = 0.0005)

B.1.1.7_B.1.177_meanrate0.0005 = 
  unpack_two_Lineages(ofn = "d1-tN-sampler1_B.1.1.7_vs_matchSample_B.1.177_meanrate5e-04.rds",
                      Lineage_main = 'B.1.1.7',
                      Lineage_matched = 'B.1.177',
                      dedup = "",
                      meanrate = 0.0005)

B.1.1.7_B.1.177_deduped_meanrate0.0005 = 
  unpack_two_Lineages(ofn = "d1-tN-sampler1_B.1.1.7_vs_matchSample_B.1.177_meanrate5e-04-dedup.rds",
                      Lineage_main = 'B.1.1.7',
                      Lineage_matched = 'B.1.177',
                      dedup = "-dedup",
                      meanrate = 0.0005)


result = do.call(rbind, lapply(list(B.1.1.7_matched_notB.1.1.7_meanrate0.001,
                                    B.1.1.7_matched_notB.1.1.7_deduped_meanrate0.001 ,
                                    B.1.1.7_B.1.177_meanrate0.001 ,
                                    B.1.1.7_B.1.177_deduped_meanrate0.001 ,
                                    B.1.1.7_matched_notB.1.1.7_meanrate0.0005 ,
                                    B.1.1.7_matched_notB.1.1.7_deduped_meanrate0.0005, 
                                    B.1.1.7_B.1.177_meanrate0.0005 ,
                                    B.1.1.7_B.1.177_deduped_meanrate0.0005 ), function(x) {
                                      
                                      rbind(
                                        c(analysis = x$ofn, outcome ='Reproduction number B.1.1.7', 
                                          round(apply( x$q_gr[ x$time >  decimal_date( as.Date('2020-11-01')) , ], 2 , function(x) median(na.omit(x)) ), 2)),
                                        
                                        c(analysis = " ", outcome ='Reproduction number control',
                                          round(apply( x$q_magr[ x$time >  decimal_date( as.Date('2020-11-01')) , ], 2 , function(x) median(na.omit(x)) ), 2)),
                                        
                                        
                                        c(analysis = " ", outcome ='Ratio reproduction numbers B.1.1.7:control',
                                          round(apply( x$qRr[ x$time >  decimal_date( as.Date('2020-11-01')) , ], 2 , function(x) median(na.omit(x)) ), 2))
                                        
                                      )
                                    }
)
)



write.csv(result, 'results.csv')

























list.files( "C:/Users/lilyl/OneDrive/Documents/N501Y_analyses/results/2021-01-18")
