#'
#' Functions required to output scanner markdown figures
#'  


#'
#'
#'dedup function

dedup_data = function(s = PHE)



#' freq_plots


freq_plots = function(sdf = csd_select
                      , kp_nodes = kp_nodes
                      , growth_rank = growth_rank2){ 
  
  library(zoo)
  library(tidyr)
  #growth_rank = readRDS(glue('{path_to_scanner}growth_rank.rds'))
  #kp_nodes = readRDS(glue('{path_to_scanner}kp_nodes.rds'))
  
  
  if(!is.null(mutations) & is.null(lineages)){ 
    
    #mutation cumulative plots 
    
    sdf_freq = data.frame(table(sdf$sel_mutation, sdf$sample_date))
    
    sdf_list = list()
    j = 0 
    for(i in mutations) { 
      j = j+1
      sdf_list[[j]] = sdf_freq%>%
        filter(as.character(Var1) == i) %>% 
        mutate(date = ymd(Var2)) %>%
        mutate(cumsum = cumsum(Freq)) %>%
        complete(date = full_seq(date, period = 1), fill = list(Freq = 0)) %>%
        mutate(cum_rolling10 = rollapplyr(Freq, width = 7, FUN = sum, partial = TRUE)) %>%
        drop_na(cumsum) 
      
    }
    
    sdf_freq = dplyr::bind_rows(sdf_list)
    # mutation plots 
    
    
    plcr = ggplot(sdf_freq, aes(date, cum_rolling10, group = Var1, colour = Var1)) + geom_line() + 
      theme_classic() + 
      labs(x = "Sample date (7 day rolling sum)", y = "Frequency of mutation (7 day cumulative sum)", color = "Mutation") + 
      scale_color_brewer(palette = "Dark2")
    
    sdf_freq$Var1 = as.character(sdf_freq$Var1)
    sdf_freq = sdf_freq[!is.na(sdf_freq$Var1),]
    
    plcs = ggplot(sdf_freq, aes(date, cumsum, group = as.character(Var1), 
                                colour = as.character(Var1))) + geom_line() + 
      theme_classic() + 
      labs(x = "Sample date", y = "Frequency of mutation (cumulative sum)", color = "Mutation") + 
      scale_color_brewer(palette = "Dark2")
    
    
    # cluster frequency plots
    sdf = csd_cluster
    sdf = csd_cluster[csd_cluster$cluster_id %in% kp_nodes,]
    
    for(i in sdf$cluster_id){ 
      sdf$growth_rank2[sdf$cluster_id == i] = growth_rank[names(growth_rank) == i]
      
    }
    
    sdf_freq = data.frame(table(sdf$growth_rank2, sdf$sample_date))
    
    sdf_list = list()
    j = 0 
    for(i in unique(sdf$growth_rank2)) { 
      j = j+1
      sdf_list[[j]] = sdf_freq%>%
        filter(as.character(Var1) == i) %>% 
        mutate(date = ymd(Var2)) %>%
        mutate(cumsum = cumsum(Freq)) %>%
        complete(date = full_seq(date, period = 1), fill = list(Freq = 0)) %>%
        mutate(cum_rolling10 = rollapplyr(Freq, width = 7, FUN = sum, partial = TRUE)) %>%
        drop_na(cumsum) 
      
    }
    
    sdf_freq = dplyr::bind_rows(sdf_list)
    
    
  }
  
  if(!is.null(lineages) & is.null(mutations)) { 
    
    sdf = csd_select
    sdf$lineage_muts = as.character(sdf$lineage_muts)
    sdf_freq = data.frame(table(sdf$lineage_muts, ymd(sdf$sample_date)))
    
    sdf_list = list()
    j = 0 
    for(i in lineages) { 
      j = j+1
      sdf_list[[j]] = sdf_freq%>%
        filter(as.character(Var1) == i) %>% 
        mutate(date = ymd(Var2)) %>%
        mutate(cumsum = cumsum(Freq)) %>%
        complete(date = full_seq(date, period = 1), fill = list(Freq = 0)) %>%
        mutate(cum_rolling10 = rollapplyr(Freq, width = 7, FUN = sum, partial = TRUE)) %>%
        drop_na(cumsum)
    }
    
    sdf_freq = dplyr::bind_rows(sdf_list)
    
    plcr = ggplot(sdf_freq, aes(date, cum_rolling10, group = Var1, colour = Var1)) + geom_line() + 
      theme_classic() + 
      labs(x = "Sample date (7 day rolling sum)", y = "Frequency of lineage (7 day cumulative sum)", color = "Lineage") + 
      scale_color_brewer(palette = "Dark2")
    
    sdf_freq$Var1 = as.character(sdf_freq$Var1)
    sdf_freq = sdf_freq[!is.na(sdf_freq$Var1),]
    
    plcs = ggplot(sdf_freq, aes(date, cumsum, group = as.character(Var1), 
                                colour = as.character(Var1))) + geom_line() + 
      theme_classic() + #geom_point() + + 
      labs(x = "Sample date", y = "Frequency of lineage (cumulative sum)", color = "Lineage") + 
      scale_color_brewer(palette = "Dark2")
    
    
    # cluster frequency plots
    sdf = csd_cluster
    sdf = csd_cluster[csd_cluster$cluster_id %in% kp_nodes,]
    
    for(i in sdf$cluster_id){ 
      sdf$growth_rank2[sdf$cluster_id == i] = growth_rank[names(growth_rank) == i]
      
    }
    
    sdf_freq = data.frame(table(sdf$growth_rank2, sdf$sample_date))
    
    sdf_list = list()
    j = 0 
    for(i in unique(sdf$growth_rank2)) { 
      j = j+1
      sdf_list[[j]] = sdf_freq%>%
        filter(as.character(Var1) == i) %>% 
        mutate(date = ymd(Var2)) %>%
        mutate(cumsum = cumsum(Freq)) %>%
        complete(date = full_seq(date, period = 1), fill = list(Freq = 0)) %>%
        mutate(cum_rolling10 = rollapplyr(Freq, width = 7, FUN = sum, partial = TRUE)) %>%
        drop_na(cumsum) 
      
    }
    
    sdf_freq = dplyr::bind_rows(sdf_list)
  }
  
  
  plcr2 = ggplot(sdf_freq, aes(date, cum_rolling10, group = Var1, colour = Var1)) + geom_line() + 
    theme_classic() + 
    labs(x = "Sample date (7 day rolling sum)", y = "Frequency of cluster (7 day cumulative sum)", color = "Mutation") + 
    scale_color_brewer(palette = "Dark2")
  
  
  plcs2 = ggplot(sdf_freq, aes(date, cumsum, group = as.character(Var1), 
                               colour = as.character(Var1))) + geom_line() + 
    theme_classic() + #geom_point() +
    #scale_x_date(limits = c(ymd("2021-01-01"), ymd("2021-01-31"))) + 
    labs(x = "Sample date", y = "Frequency of cluster (cumulative sum)", color = "Lineage, Growth rank") + 
    scale_color_brewer(palette = "Dark2")
  
  
  ggsave(filename = glue('{path_to_scanner}cum_rolling_{max_date}_{min_date}.png'), plcr, width = 6.5, height = 3.5)
  ggsave(filename = glue('{path_to_scanner}cum_samples_{max_date}_{min_date}.png'), plcs, width = 6.5, height = 3.5)
  ggsave(filename = glue('{path_to_scanner}cum_rolling_clusters_{max_date}_{min_date}.png'), plcr2, width = 6.5, height = 3.5)
  ggsave(filename = glue('{path_to_scanner}cum_samples_clusters_{max_date}_{min_date}.png'), plcs2, width = 6.5, height = 3.5)
  
}


#'
#'GAM outputs 
#'

variable_frequency_day_report <- function(s
                                          , variable='type'
                                          , value='clade'
                                          ,  mint = -Inf 
                                          , maxt = Inf
                                          , detailed=TRUE 
                                          , form = y~s(sample_time, bs = 'gp', k = 5) )
{
  library( lubridate ) 
  stopifnot( is.numeric( s$sample_time ) )
  s <- s[ s$sample_time >= mint  &  s$sample_time <= maxt , ]
  rownames(s) <- NULL
  s <- s[ order(s$sample_time ) , ]
  
  if ( 'weight' %in% colnames(s)){
    s$weight <- s$weight / mean( s$weight ) 
  } else{
    s$weight <- 1
  }
  
  ss <- split( s, s[[variable]]==value )
  
  library( mgcv ) 
  .s1 <- s[ order( s$sample_time ) , ]
  .s1 <- .s1[ .s1$sample_time >= (min(ss[['TRUE']]$sample_time)-1/365), ]
  .s1$y <- .s1[[variable]]==value
  m = mgcv::gam( form , family = binomial(link='logit') , data = .s1 )
  .s1$estimated = predict( m )
  estdf = data.frame( time = .s1$sample_time, estimated_logodds = .s1$estimated )
  estdf <- estdf[ !duplicated(estdf$time), ]
  estdf <- estdf[ order( estdf$time ) , ]
  estdf <- cbind( estdf, t( sapply( 1:nrow(estdf), function(k){
    .s2 <- .s1[ .s1$sample_time == estdf$time[k] , ]
    n <- nrow( .s2 )
    lo= estdf$estimated_logodds[k] 
    f = exp( lo) / ( 1 + exp( lo ))
    ub = qbinom( .975, size = n, prob = f  )  / n
    lb = qbinom( .025, size = n, prob = f  )  /n
    n1 <- sum( .s2[[variable]]==value)
    n2 <- n - n1 
    ef = n1 / n 
    c(lb= log(lb/(1-lb)), ub = log(ub/(1-ub)) , n = n, weights = f*(1-f)/n, logodds = log(ef / ( 1 - ef )) )
  })) )
  
  
  estdf
  
}




logistic_growth_stat2 <- function(u, e1 = envs, generation_time_scale = 6.5/365)
{
  ta = get_comparator_sample(u) 
  if ( is.null( ta ))
    return(0)
  
  tu = e1$descendantSids[[u]] 
  sta = na.omit(e1$sts[ na.omit(ta) ])
  stu = e1$sts[ na.omit( e1$descendantSids[[u]] )  ]
  ta = ta[ta %in% names(sta)]
  tu = tu[tu %in% names(stu)]
  X = data.frame( sample_time = c( sta, stu ), type = c( rep('control', length(ta)), rep('clade',length(tu)) ), 
                  sequence_name = c(ta,tu) )
  X$weights = ifelse(X$sequence_name %in% amd$sequence_name,
                     amd$coverage_weight[amd$sequence_name %in% c(ta,tu)], NA )
  X$lineage = growth_rank[names(growth_rank) == u]
  X
  
}



gam_plots = function(e1 = envs
                     , kp_nodes = kp_nodes
                     , growth_rank = growth_rank2) {
  #growth_rank = readRDS(glue('{path_to_scanner}growth_rank.rds'))
  #kp_nodes = readRDS(glue('{path_to_scanner}kp_nodes.rds'))
  
  
  
  us = as.integer(kp_nodes)
  Xs = lapply(us, logistic_growth_stat2)
  ss = lapply(Xs, variable_frequency_day_report)
  names(ss) = us
  
  for(i in us){ 
    
    ss[[as.character(i)]]$grouping = growth_rank[names(growth_rank) == i]
    
  }
  
  estdf = dplyr::bind_rows(ss)
  
  saveRDS(estdf, glue('{path_to_scanner}_estdf_{max_date}_{min_date}.rds'))
  #, weight=weights
  
  pl = ggplot(data= estdf[estdf$time > 2020.88,]) + 
    geom_point( aes( x = as.Date( date_decimal( time ) ), y = logodds, color = grouping, size=weights ))  + 
    theme_minimal() + facet_wrap(~lineage, scales = "free_x") + scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    theme( legend.pos='' ) +
    xlab('') + ylab("Log odds sample") +
    geom_path( aes(x=as.Date(date_decimal( time )), y=estimated_logodds,color = grouping ), size = 1) + 
    geom_ribbon( aes(x = as.Date(date_decimal( time )), ymin = lb, ymax = ub, fill = grouping),   alpha = .25 )
  
  
  ggsave(filename = glue('{path_to_scanner}logodds_{max_date}_{min_date}.png'), pl, width = 6.5, height = 3.5)
  
}  



#'
#'
#'
#'plot_maps

plot_maps = function(sdf = csd_cluster
                     , path_to_data = path_to_data){
  
  path_to_map = glue("{path_to_data}UK_NUTS_1_boundaries.rds")
  
  gg = readRDS(path_to_map)
  sdf = csd_cluster
  sdf = sdf %>% filter(cluster_id %in% kp_nodes)
  
  
  #grouping label based on scanner selection
  
  if(!is.null(lineages) & is.null(mutations)){ 
    labs_color = "Lineage, Growth rank"
    sdf$growth_rank2 = sdf$growth_rank
  }
  if(!is.null(mutations)){ 
    labs_color = "Lineage, Growth rank"}
  
  if(is.null(lineages) & is.null(mutations)){ 
    labs_color = "Growth rank"}
  
  sdf = sdf %>% select(longitude, latitude, logistic_growth_rate, growth_rank, growth_rank2, sample_date) 
  sdf$logistic_growth_rate = as.numeric(sprintf(as.numeric(sdf$logistic_growth_rate), fmt = '%#.2f')) 
  gg1 = gg + geom_point(data = sdf, aes(x = longitude, y = latitude, group = logistic_growth_rate
                                        , color= logistic_growth_rate), size = 2, alpha = 0.5) +
    scale_color_gradient2(midpoint = mean(sdf$logistic_growth_rate)
                          , low = "#00AFBB", mid = "#E7B800",high = "#FC4E07", space = "Lab" ) +
    labs(color = "Logistic growth rate")
  
  gg2 = gg + geom_point(data = sdf, aes(x = longitude, y = latitude, group = growth_rank2, color = growth_rank2)
                        ,size = 2, alpha = 0.5) +
    scale_color_brewer(palette = "Dark2") + labs(color = labs_color)
  
  
  gg3 = gg + geom_point(data = sdf, aes(x = longitude, y = latitude, color= sample_date)
                        ,size = 2, alpha = 0.4)  +
    labs(color = "Sample date")
  
  ggsave(glue('{path_to_scanner}log_growth_map_{max_date}_{min_date}.png'), gg1, width = 6.5, height = 3.5)
  ggsave(glue('{path_to_scanner}region_map_{max_date}_{min_date}.png'), gg2, width = 6.5, height = 3.5)
  ggsave(glue('{path_to_scanner}sample_time_map_{max_date}_{min_date}.png'), gg3, width = 6.5, height = 3.5)
  
}



#'
#'
#'
#'
#'mutation_plots




#' @export 
cluster_muts_new = function( scanner_env 
                             , nodes 
                             , mutfn = list.files( paste0( path_to_data) ,   patt = 'cog_global_[0-9\\-]+_mutations.csv' , full.name = TRUE)
                             , min_seq_contrast = 1 # sequences in ancestor clade
                             , overlap_threshold = .9
)
{
  e1 = as.environment( scanner_env )
  attach( e1 )
  
  Ygr1 = Y[ Y$node_number %in% nodes , ] 
  
  mdf = read.csv( mutfn, stringsAs=FALSE )
  ku = 1
  cms = lapply( 1:nrow(Ygr1) , function(ku) {
    mdf.u = mdf[ mdf$sequence_name %in% strsplit( Ygr1$tips[ku], split = '\\|')[[1]] , ]
    u = nodes[ku] 
    
    vtabu = sort( table( do.call( c, strsplit( mdf.u$variants, split='\\|' )  ) ) / nrow( mdf.u ) )
    
  })
  
  detach(e1)
  cms
  
}


#'   mutation plots 
#'   
#'   

mutations_plots = function(cmts2 = cmts2
                           , kp_nodes = kp_nodes
                           , prop_mutations = prop_mutations) {
  
  cms_data = list()
  for(i in kp_nodes){
    
    cms_data[[as.character(i)]] =  data.frame(  mutants = c(names(cmts2[[as.character(i)]])),
                                                props = c(cmts2[[as.character(i)]]),
                                                cluster_id = c(rep(i,length(cmts2[[as.character(i)]]))), 
                                                gr_muts = c(rep(growth_rank[names(growth_rank) == i], length = length(cmts2[[as.character(i)]])))
    )
    
  }
  
  
  
  
  
  cms_data = bind_rows(cms_data)
  
  
  cms_data_syn = cms_data %>% 
    filter(str_detect(mutants, "synSNP"))
  
  cms_data_s = cms_data %>% 
    filter(str_detect(mutants, "S:"))
  
  cms_data_n = cms_data %>% 
    filter(str_detect(mutants, "N:"))
  syns_out = cms_data_syn$mutants
  
  
  plot <- ggplot(cms_data[cms_data$mutants %nin%  syns_out & cms_data$props >= prop_mutations,], aes(x=reorder(mutants, props), props, fill=as.factor(gr_muts)))
  plot <- plot + geom_bar(stat = "identity", position = 'dodge') + theme_classic()+ facet_wrap(~gr_muts, scale = "free") + scale_fill_brewer(palette = "Dark2")
  plot <- plot + labs(y = "Proportion of sequences", x = "Mutations", fill = "Lineage") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  plot <- plot + theme(legend.position = "none")
  
  ggsave(filename = glue('{path_to_scanner}mutations_{max_date}_{min_date}.png'), plot, width = 6, height = 7.5 )
  
}

#~ prop_mutations = 0.25
#~ mutations_plots()

clock_outlier_plot = function(){ 
  
  for(i in kp_nodes){ 
    
    Ygr2$growth_rank[Ygr2$cluster_id == as.character(i)] = growth_rank[names(growth_rank) == i]
    
  }
  
  cplot = ggplot(Ygr2, aes(logistic_growth_rate, clock_outlier, color = growth_rank, size = cluster_size)) + geom_point() + theme_minimal() + 
    scale_color_brewer(palette = "Dark2") + labs(y = "Clock outlier", x = "Logistic growth rate", color = "Lineage, Mutant, Log rank", size = "Cluster size")
  
  ggsave(filename = glue('{path_to_scanner}clockoutlier_{max_date}_{min_date}.png'), plot, width = 7.5, height = 4)  
  
  
}

