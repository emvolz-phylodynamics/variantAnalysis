#'Table for scanner report 
#'
#'
#'
#'@requires lineages, mutations, cut_off, sequence scanner data, Ygr1

report_table = function(csd = s, lineages, mutations, cut_off, Ygr1) { 

  
  #just lineages 
  
  if(!is.null(lineages) & is.null(mutations)){  
    
    csd_lin = list()
    csd_names = NA
    csd_all  = list() 
    
    for(i in 1:length(lineages)){ 
      csd_all[[i]] = csd[csd$lineage_muts %in% lineages[[i]],]
      csd_all[[i]]$sel_lineage = lineages[[i]]
      csd_all[[i]]$sel_log_rank = glue('{csd_all[[i]]$sel_lineage},{csd_all[[i]]$log_rank}')
    }
    
    csd_all1 = csd_all
    names(csd_all1) = lineages
    csd_all = dplyr::bind_rows(csd_all)
    
    for(i in unique(csd_all$log_rank)){ 
      
      csd_all$sel_lineage[csd_all$log_rank == i] = paste(unique(csd_all$sel_lineage[csd_all$log_rank == i]), collapse = ' ')
      
    }
    csd_all$growth_rank = glue('{csd_all$sel_lineage}, {csd_all$log_rank}')
    growth_rank = unique(csd_all$growth_rank)
    gr_names=NA
    
    csd_select = csd_all
    
    for(i in 1:length(growth_rank)){
      j = as.integer(sub(".*, ", "", growth_rank[i]))
      gr_names[i]  = unique(as.character(csd_all$cluster_id[csd_all$log_rank == j]))
    }
    names(growth_rank) = gr_names
    growth_rank2 = growth_rank
    sel_log_rank = unique(csd_all$sel_log_rank)
    log_rank = sub(".*,","", sel_log_rank)
    sel_lineage = sub(",.*", "", sel_log_rank) 
    
    for(j in 1:length(log_rank)){
      csd_lin[[j]] = csd[csd$log_rank %in% log_rank[[j]],]
      csd_lin[[j]]$sel_lineage = sel_lineage[[j]]
      csd_lin[[j]]$sel_lineage_pres = ifelse(csd_lin[[j]]$sequence_name %in% csd_all1[[sel_lineage[j]]]$sequence_name, 
                                             TRUE, FALSE)
      
      csd_names[j] = log_rank[j]
      
    }
    
    names(csd_lin) = csd_names
    
    csd_lin_all = dplyr::bind_rows(csd_lin)
    for(i in csd_lin_all$cluster_id){
      csd_lin_all$growth_rank[csd_lin_all$cluster_id == i] = growth_rank[names(growth_rank) == i]
    }
    csd_cluster = csd_lin_all
    
    perc_all = list() 
    #cluster_ids = unique(csd$log_rank[ csd$lineage_muts %in% lineages])
    for(j in 1:length(lineages)){ 
      perc_lin = list() 
      cluster_ids = unique(csd_all1[[j]]$log_rank)
      
      for(i in 1:length(cluster_ids)) { 
        
        perc_lin[[i]] = data.frame(table(csd_lin_all$sel_lineage_pres[csd_lin_all$sel_lineage == lineages[j] & 
                                                                        csd_lin_all$log_rank == cluster_ids[i]]))
        perc_lin[[i]]$percent = perc_lin[[i]]$Freq/sum(perc_lin[[i]]$Freq) * 100
        perc_lin[[i]]$log_rank = cluster_ids[i]
        perc_lin[[i]]$sel_lineages = lineages[[j]]
        perc_lin[[i]] = perc_lin[[i]][perc_lin[[i]]$Var1 == TRUE,]
        
      }
      perc_lin = dplyr::bind_rows(perc_lin)
      perc_all[[j]] = perc_lin
      
    }
    
    names(perc_all) = lineages
    cluster_ids = log_rank
    
    #defining mutations 
    Ygr2 = Ygr1[Ygr1$log_rank %in% cluster_ids,]
    kp_nodes = as.numeric(unique(Ygr2$cluster_id))
    max_growth_rate = sprintf(as.numeric(max(Ygr2$logistic_growth_rate)), fmt = '%#.2f') 
    n_clusters_included = length(cluster_ids)
    cmts1 = cluster_muts(scanner_env = envs, nodes = kp_nodes, overlap_threshold = defining_mutations_cut_off/100)
    cmts2 = cluster_muts_all(scanner_env = envs, nodes = kp_nodes, overlap_threshold = defining_mutations_cut_off/100)
    cdels1 = cluster_dels(scanner_env = envs, nodes = kp_nodes, overlap_threshold = defining_mutations_cut_off/100)
    cdels2 = cluster_dels_all(scanner_env = envs, nodes = kp_nodes, overlap_threshold = defining_mutations_cut_off/100)
    names(cmts1) = sort(kp_nodes)
    names(cmts2) = sort(kp_nodes)
    #remove synsnps 
    muts = list()
    muts_all = list()
    for(i in 1:length(Ygr2$log_rank)) { 
      
      cmts1[[i]] = cmts1[[i]][!str_detect(cmts1[[i]],pattern="synSNP")]
      muts[[i]] = cmts1[[i]][str_detect(cmts1[[i]],pattern= moc)]

      cmts3 = names(cmts2[[i]])
      muts_all[[i]] = cmts3[str_detect(cmts3,pattern= moc)]

    }
    
    names(muts) = kp_nodes
    names(muts_all) = kp_nodes
    
    Ygr2 = list() 
    
    for(j in 1:length(lineages)){ 
      
      Ygr2[[j]] = Ygr1[Ygr1$log_rank %in% csd_all1[[j]]$log_rank,]
      cluster_ids = unique(csd_all1[[j]]$log_rank)
      
      for(i in cluster_ids){ 

        Ygr2[[j]]$percent[Ygr2[[j]]$log_rank == i] = perc_all[[j]]$percent[perc_all[[j]]$log_rank == i]
        Ygr2[[j]]$sel_lineage[Ygr2[[j]]$log_rank == i] = perc_all[[j]]$sel_lineages[perc_all[[j]]$log_rank == i]
        Ygr2[[j]]$tips = NA
        
        if(Ygr2[[j]]$percent[Ygr2[[j]]$log_rank == i] < 100){ 
          
          other_perc = glue('({as.integer(table(csd_lin[[as.character(i)]]$lineage_muts[csd_lin[[as.character(i)]]$lineage_muts 
                                                                     %nin% unique(Ygr2[[j]]$sel_lineage)])/nrow(csd_lin[[as.character(i)]])*100)})')
          other_lin = names(table(csd_lin[[as.character(i)]]$lineage_muts[csd_lin[[as.character(i)]]$lineage_muts 
                                                                          %nin% unique(Ygr2[[j]]$sel_lineage)])/nrow(csd_lin[[as.character(i)]]))
          
          Ygr2[[j]]$other_lineages[Ygr2[[j]]$log_rank == i] = 
            paste(glue('{other_lin}{other_perc}'), 
                  collapse = ' ')
        } else { 
          
          Ygr2[[j]]$other_lineages[Ygr2[[j]]$log_rank == i] = ""
          
        }
        
      }
    }
    
    
  
    if(!is.null(cut_off)){ 
      
      if(cut_off < 5){
        
        for(i in 1:length(lineages)) { 
          
          Ygr2[[i]] = Ygr2[[i]][Ygr2[[i]]$logistic_growth_rate >= cut_off,]
          
        } 
        
      } else { 
        
        for(i in 1:length(lineages)){ 
          
          Ygr2[[i]] = Ygr2[[i]][order(Ygr2[[1]]$log_rank),]
          Ygr2[[i]] = Ygr2[[i]][1:cut_off,]
          
        }
        
        
      }
    }
    
    Ygr2 = dplyr::bind_rows(Ygr2)
    
    for(i in kp_nodes){ 
      
      Ygr2$defining_muts[Ygr2$cluster_id == i] = paste(cmts1[[as.character(i)]], collapse = ' ') 
      Ygr2$defining_dels[Ygr2$cluster_id == i] = paste(cdels1[[as.character(i)]], collapse = ' ')
      Ygr2$muts_in[Ygr2$cluster_id == i] = ifelse(is_empty(muts[[as.character(i)]]) & 
                                                  is_empty(muts_all[[as.character(i)]]), glue('{moc} Absent'), 
                                                ifelse(is_empty(muts[[as.character(i)]]) & muts_all[[as.character(i)]] == moc, glue('{moc} Present'), glue('{moc} Defining')))
      
    }
    
    Ygr2 = Ygr2[Ygr2$logistic_growth_rate > 0,]
    Ygr2$logistic_growth_rate = sprintf(as.numeric(Ygr2$logistic_growth_rate), fmt = '%#.2f')  
    Ygr2$percent = as.integer(Ygr2$percent) 
    Ygr2 = Ygr2[Ygr2$percent > prop_cluster,]
    Ygr2$cluster_size_perc = glue('{Ygr2$cluster_size}({Ygr2$percent})')
    Ygr2$grouping = glue('{Ygr2$sel_lineage},{Ygr2$log_rank}')
    
    
    Ygr2_output = Ygr2 %>% select(sel_lineage, cluster_size_perc, other_lineages,  most_recent_tip, least_recent_tip, logistic_growth_rate, log_rank, muts_in, defining_dels, defining_muts)
    
    colnames(Ygr2_output) = c("Selection", "Size (%)", "Lineages (%)","Most recent tip", 
                              "Least recent tip", "Log growth rate", "Growth rate rank", 
                              "Mutation status", "Defining deletions", "Defining mutations")
    rownames(Ygr2_output) <- c()

    
  }
  
  
  # just mutations 
  if(is.null(lineages) & !is.null(mutations)){  
    
    csd_muts = list()
    
    for(i in 1:length(mutations)){ 
      csd_muts[[i]] = csd %>% filter(str_detect(variants, mutations[[i]]))
      csd_muts[[i]]$sel_mutation = mutations[[i]]
      csd_muts[[i]]$sel_mut_log_rank = glue('{csd_muts[[i]]$sel_mutation},{csd_muts[[i]]$log_rank}')

    }
  
    
    csd_muts1 = csd_muts
    names(csd_muts1) = mutations
    csd_muts = dplyr::bind_rows(csd_muts)
   
    
    for(i in unique(csd_muts$log_rank)){ 
      
      csd_muts$sel_muts_all[csd_muts$log_rank == i] = paste(unique(csd_muts$sel_mutation[csd_muts$log_rank == i]), collapse = ' ')
      csd_muts$sel_lineage[csd_muts$log_rank == i] = paste(unique(csd_muts$lineage_muts[csd_muts$log_rank == i]), collapse = ' ')
      
    }
    
    csd_muts$growth_rank = glue('{csd_muts$sel_lineage}, {csd_muts$sel_muts_all}, {csd_muts$log_rank}')
    csd_muts$growth_rank2 = glue('{csd_muts$sel_lineage}, {csd_muts$log_rank}')
    csd_select = csd_muts
    
    growth_rank = unique(csd_muts$growth_rank)
    growth_rank2 = unique(csd_muts$growth_rank2)
    
    gr_names=NA
    for(i in 1:length(growth_rank)){
    j = as.integer(sub(".*, ", "", growth_rank[i]))
    gr_names[i]  = unique(as.character(csd_muts$cluster_id[csd_muts$log_rank == j]))
    }
    names(growth_rank) = gr_names
    names(growth_rank2) = gr_names
    
    sel_mut_log_rank =unique(csd_muts$sel_mut_log_rank)
    log_rank = sub(".*,","", sel_mut_log_rank)
    sel_mutations = sub(",.*", "", sel_mut_log_rank) 
    
    csd_lin = list()
    csd_names = NA
    
    for( j in 1:length(log_rank)) {
      csd_lin[[j]] = csd[csd$log_rank %in% log_rank[[j]],] 
      csd_lin[[j]]$sel_mutation = sel_mutations[[j]]
      csd_lin[[j]]$sel_mutation_pres = ifelse(csd_lin[[j]]$sequence_name %in% csd_muts1[[sel_mutations[j]]]$sequence_name, 
                                              TRUE, FALSE)
      csd_names[j] = log_rank[[j]]
      
      
    }
    
    names(csd_lin) = csd_names
    
    csd_lin_all = dplyr::bind_rows(csd_lin)
    
    for(i in unique(csd_lin_all$cluster_id)){
      csd_lin_all$growth_rank[csd_lin_all$cluster_id == i] = growth_rank[names(growth_rank) == i]
      csd_lin_all$growth_rank2[csd_lin_all$cluster_id == i] = growth_rank2[names(growth_rank2) == i]
    }
    
    csd_cluster = csd_lin_all
  
    perc_all = list() 
    
    for(j in 1:length(mutations)) { 
      
      perc_lin = list() 
      cluster_ids = unique(csd_muts1[[j]]$log_rank)
      
      for(i in 1:length(cluster_ids)) { 
        
        perc_lin[[i]] = data.frame(table(csd_lin_all$sel_mutation_pres[csd_lin_all$sel_mutation == mutations[j] & 
                                                                         csd_lin_all$log_rank == cluster_ids[i]]))
        perc_lin[[i]]$percent = perc_lin[[i]]$Freq/sum(perc_lin[[i]]$Freq) * 100
        perc_lin[[i]]$log_rank = cluster_ids[i]
        perc_lin[[i]]$sel_mutations = mutations[[j]]
        perc_lin[[i]] = perc_lin[[i]][perc_lin[[i]]$Var1 == TRUE,]
        
      }
      perc_lin = dplyr::bind_rows(perc_lin)
      perc_all[[j]] = perc_lin
      
    } 
    
    names(perc_all) = mutations
    cluster_ids = log_rank
    
    #defining mutations 
    Ygr2 = Ygr1[Ygr1$log_rank %in% cluster_ids,]
    kp_nodes = as.numeric(unique(Ygr2$cluster_id))
    max_growth_rate = sprintf(as.numeric(max(Ygr2$logistic_growth_rate)), fmt = '%#.2f') 
    n_clusters_included = length(cluster_ids)
    cmts1 = cluster_muts(scanner_env = envs, nodes = kp_nodes, overlap_threshold = defining_mutations_cut_off/100)
    cmts2 = cluster_muts_all(scanner_env = envs, nodes = kp_nodes, overlap_threshold = defining_mutations_cut_off/100)
    cdels1 = cluster_dels(scanner_env = envs, nodes = kp_nodes, overlap_threshold = defining_mutations_cut_off/100)
    cdels2 = cluster_dels_all(scanner_env = envs, nodes = kp_nodes, overlap_threshold = defining_mutations_cut_off/100)
    names(cmts1) = sort(kp_nodes)
    names(cmts2) = sort(kp_nodes)
    names(cdels1) = sort(kp_nodes)
    names(cdels2) = sort(kp_nodes)

    for(i in 1:length(Ygr2$log_rank)) { 
      
      cmts1[[i]] = cmts1[[i]][!str_detect(cmts1[[i]],pattern="synSNP")]
      cmts3 = names(cmts2[[i]])
    }
    
    dpa = list()
    dpa2 = list() 
    dpa3 = list()
    dpa4 = list() 
    for(i in 1:length(mutations)){
      
      for(j in 1:length(kp_nodes)){
        
        cmts3 = names(cmts2[[j]])
        dpa[[j]] = cmts3[str_detect(cmts3,pattern = mutations[[i]])] 
        dpa2[[j]] = cmts1[[j]][str_detect(cmts1[[j]],pattern = mutations[[i]])]
        
      }
      names(dpa) = sort(kp_nodes)
      names(dpa2) = sort(kp_nodes)
      dpa3[[i]] = dpa # present 
      dpa4[[i]] = dpa2 # defining
    }
    
    #setting up defining mutants, lineages, other lineages
    Ygr2 = list()
    
    for(j in 1:length(mutations)){
      
      Ygr2[[j]] = Ygr1[Ygr1$log_rank %in% csd_muts1[[j]]$log_rank,] 
      cluster_ids = unique(csd_muts1[[j]]$log_rank)
      
      for(i in cluster_ids) {
        
        Ygr2[[j]]$percent[Ygr2[[j]]$log_rank == i] = perc_all[[j]]$percent[perc_all[[j]]$log_rank == i]
        lins_perc = glue('({as.integer(table(csd_lin[[as.character(i)]]$lineage_muts)/nrow(csd_lin[[as.character(i)]])*100)})')
        lins = names(table(csd_lin[[as.character(i)]]$lineage_muts))
        Ygr2[[j]]$lineages[Ygr2[[j]]$log_rank == i] = paste(glue('{lins}{lins_perc}'), collapse = ' ')
        Ygr2[[j]]$sel_mutation[Ygr2[[j]]$log_rank == i] = perc_all[[j]]$sel_mutations[perc_all[[j]]$log_rank == i]
        Ygr2[[j]]$tips = NA
      }
    }
    
    for(i in sort(kp_nodes)) { 
      
      for(j in 1:length(mutations)) { 
        
        Ygr2[[j]]$defining_muts[Ygr2[[j]]$cluster_id == i] = paste(cmts1[[as.character(i)]], collapse = ' ') 
        Ygr2[[j]]$defining_dels[Ygr2[[j]]$cluster_id == i] = paste(cdels1[[as.character(i)]], collapse = ' ') 
        
        Ygr2[[j]]$dpa[Ygr2[[j]]$cluster_id == i] = ifelse(!is_empty(dpa4[[j]][[as.character(i)]]) & 
                                                            !is_empty(dpa3[[j]][[as.character(i)]]), "Defining", 
                                                          ifelse(is_empty(dpa4[[j]][[as.character(i)]]) &
                                                                   !is_empty(dpa3[[j]][[as.character(i)]]) , "Present", 
                                                                 "Absent"))
      }
      
    }
  
    if(!is.null(cut_off)){ 
      
      if(cut_off < 5){
        
        for(i in 1:length(mutations)) { 
          
          Ygr2[[i]] = Ygr2[[i]][Ygr2[[i]]$logistic_growth_rate >= cut_off,]
          
        } 
        
      } else { 
        
        for(i in 1:length(mutations)){ 
          
          Ygr2[[i]] = Ygr2[[i]][Ygr2[[i]]$log_rank <= cut_off,]
          
        }
        
        
      }
    }
    
    
    Ygr2 = dplyr::bind_rows(Ygr2)
    Ygr2$logistic_growth_rate = sprintf(as.numeric(Ygr2$logistic_growth_rate), fmt = '%#.2f')
    Ygr2 = Ygr2[as.numeric(Ygr2$logistic_growth_rate) > log_growth_rate_cut_off,]
    Ygr2 = Ygr2[Ygr2$percent > prop_cluster,]
    Ygr2$percent = as.integer(Ygr2$percent)
    Ygr2$cluster_size_perc = glue('{Ygr2$cluster_size}({Ygr2$percent})')
    Ygr2$grouping = glue('{Ygr2$lineages}, {Ygr2$sel_mutation}, {Ygr2$log_rank}')
    Ygr2_output = Ygr2 %>% select(sel_mutation, cluster_size_perc, lineages,  most_recent_tip, least_recent_tip, logistic_growth_rate, log_rank, dpa, defining_muts, defining_dels)
    
    colnames(Ygr2_output) = c("Selection", "Size(%)", "Lineages(%)", "Most recent tip", "Least recent tip", "Log growth rate", "Growth rate rank", "Mutation status", "Defining mutations", "Defining deletions")
    rownames(Ygr2_output) <- c()
    
    
    
  }
  


  if(is.null(lineages) & is.null(mutations) & !is.null(cut_off)){
    
    
    if(cut_off < 5){ 
      Ygr2 = Ygr1[Ygr1$log_rank >= cut_off,]
      
    }
    
    if(cut_off > 5){ 
      Ygr2 = Ygr1[Ygr1$log_rank[1:cut_off],]
      
    }
    kp_nodes = as.numeric(unique(Ygr2$cluster_id))
    max_growth_rate = sprintf(as.numeric(max(Ygr1$logistic_growth_rate)), fmt = '%#.2f')
    max_growth_rate2 = sprintf(as.numeric(max(Ygr2$logistic_growth_rate)), fmt = '%#.2f')
    min_growth_rate2 = sprintf(as.numeric(min(Ygr2$logistic_growth_rate)), fmt = '%#.2f')
    n_clusters_all = length(Ygr1$logistic_growth_rate)
    cluster_ids = Ygr2$log_rank 
    n_clusters_included = length(cluster_ids)
    cmts1 = cluster_muts(scanner_env = envs, nodes = kp_nodes, overlap_threshold = defining_mutations_cut_off/100)
    cmts2 = cluster_muts_new(scanner_env = envs, nodes = kp_nodes, overlap_threshold = defining_mutations_cut_off/100)
    names(cmts1) = sort(kp_nodes)
    names(cmts2) = sort(kp_nodes)
    #remove synsnps 
    e484k = list()
    s494p = list()
    l18f = list()
    e484k_all = list()
    s494p_all = list()
    sl18f_all = list()
    for(i in 1:length(Ygr2$log_rank)) { 
      
      cmts1[[i]] = cmts1[[i]][!str_detect(cmts1[[i]],pattern="synSNP")]
      e484k[[i]] = cmts1[[i]][str_detect((cmts1)[[i]],pattern="E484K")]
      #l18f[[i]] = cmts1[[i]][str_detect((cmts1)[[i]],pattern="L18F")]
      #s494p[[i]] = cmts1[[i]][str_detect((cmts1)[[i]],pattern="S494P")]
      cmts3 = names(cmts2[[i]])
      e484k_all[[i]] = cmts3[str_detect(cmts3,pattern="E484K")]
      #sl18f_all[[i]] = cmts3[str_detect(cmts3,pattern="L18F")]
      #s494p_all[[i]] = cmts3[str_detect(cmts3,pattern="S494P")]
    }
    names(e484k) = sort(kp_nodes)
    names(e484k_all) = sort(kp_nodes)
    
    for(i in kp_nodes){ 
      Ygr2$defining_muts[Ygr2$cluster_id == i] = paste(cmts1[[as.character(i)]], collapse = ' ') 
      Ygr2$e484K[Ygr2$cluster_id == i] = ifelse(is_empty(e484k[[as.character(i)]]) & 
                                                  is_empty(e484k_all[[as.character(i)]]), "Absent", 
                                                ifelse(is_empty(e484k[[as.character(i)]]) & e484k_all[[as.character(i)]] =="S:E484K", "Present", "Defining"))
      
    }
    
    csd_all = list()
    #library(Hmisc)
    for(i in kp_nodes){ 
      
      csd_lin = csd[csd$cluster_id %in% i,]
      cc = table(as.character(csd_lin$lineage_muts)) # change to lineage muts?
      
      for(j in names(cc)){
        csd_lin$lin_count[as.character(csd_lin$lineage_muts) == j] = cc[names(cc) == j] # change to lineage muts?
        
        
      }
      csd_lin$perc = csd_lin$lin_count/nrow(csd_lin) * 100
      csd_lin$lineages = names(cc[rank(cc) == length(cc)])
      csd_lin$other_lineages = ifelse(length(cc) > 1, paste(names(cc[rank(cc) %nin% length(cc)]), collapse = ' '), "")
      
      csd_all[[as.character(i)]] = csd_lin
      
    }
    
    for(i in kp_nodes) {
      Ygr2$lineages[Ygr2$cluster_id == i] = unique(csd_all[[as.character(i)]]$lineages)
      Ygr2$other_lineages[Ygr2$cluster_id == i] = unique(csd_all[[as.character(i)]]$other_lineages)
      Ygr2$percent[Ygr2$cluster_id == i] = unique(csd_all[[as.character(i)]]$perc[csd_all[[as.character(i)]]$lineage_muts == 
                                                                                    Ygr2$lineages[Ygr2$cluster_id == i]]) 
      
    }
    
    
    growth_rank = glue('{Ygr2$lineages} {Ygr2$other_lineages}, {Ygr2$log_rank}')
    growth_rank = unique(csd_all$growth_rank)
    gr_names=NA
    for(i in 1:length(growth_rank)){
      j = as.integer(sub(".*, ", "", growth_rank[i]))
      gr_names[i]  = as.character(Ygr2$cluster_id[Ygr2$log_rank == j])
    }
    names(growth_rank) = gr_names
    
    
    csd_lin_all = bind_rows(csd_all)
    for(i in csd_lin_all$cluster_id){
      csd_lin_all$growth_rank[csd_lin_all$cluster_id == i] = growth_rank[names(growth_rank) == i]
    }
    
    csd_cluster = csd_lin_all
    
    Ygr2$e484k_in = Ygr2$e484K
    Ygr2$logistic_growth_rate = sprintf(as.numeric(Ygr2$logistic_growth_rate), fmt = '%#.2f')   
    Ygr2$percent = as.integer(Ygr2$percent) 
    Ygr2$cluster_size_perc = glue('{Ygr2$cluster_size}({Ygr2$percent})')
    Ygr2_output = Ygr2 %>% select(lineages, other_lineages, cluster_size_perc, most_recent_tip, least_recent_tip, logistic_growth_rate, log_rank, e484k_in, defining_muts)
    
    colnames(Ygr2_output) = c("Lineages", "Other lineages", "Cluster size (%)", "Most recent tip", "Least recent tip", "Log growth rate", "Growth rate rank", "E484K status", "Defining mutations")
    rownames(Ygr2_output) <- c()
    csd_select = NA
    
  }
  
  
 if(!is.null(lineages) & !is.null(mutations)){  
    
    csd_muts = list()
    
    for(i in 1:length(mutations)){ 
      csd_muts[[i]] = csd %>% filter(str_detect(variants, mutations[[i]]))
      csd_muts[[i]]$sel_mutation = mutations[[i]]
      csd_muts[[i]]$sel_mut_log_rank = glue('{csd_muts[[i]]$sel_mutation},{csd_muts[[i]]$log_rank}')
      csd_muts[[i]] = csd_muts[[i]][csd_muts[[i]]$lineage_muts %in% lineages,]
      csd_muts[[i]]$sel_lineages_mut_log_rank = glue('{csd_muts[[i]]$lineage_muts},{csd_muts[[i]]$sel_mutation}')
      csd_muts[[i]]$sel_lineages_mut = glue('{csd_muts[[i]]$sel_mutation},{csd_muts[[i]]$lineage_muts},{csd_muts[[i]]$log_rank}')
      csd_muts[[i]]$sel_mut_lineages_log_rank = glue('{csd_muts[[i]]$lineage_muts},{csd_muts[[i]]$sel_mutation},{csd_muts[[i]]$log_rank}')
      
      
    }
    csd_muts1 = csd_muts
    names(csd_muts1) = mutations
    csd_muts = dplyr::bind_rows(csd_muts)
    
    for(i in unique(csd_muts$log_rank)){ 
      
      csd_muts$sel_muts_all[csd_muts$log_rank == i] = paste(unique(csd_muts$sel_mutation[csd_muts$log_rank == i]), collapse = ' ')
      csd_muts$sel_lineage[csd_muts$log_rank == i] = paste(unique(csd_muts$lineage_muts[csd_muts$log_rank == i]), collapse = ' ')
      
    }
    csd_muts$growth_rank = glue('{csd_muts$sel_lineage},{csd_muts$sel_muts_all},{csd_muts$log_rank}')

    growth_rank = unique(csd_muts$growth_rank)

    gr_names=NA
    for(i in 1:length(growth_rank)){
      j = as.integer(sub(".*,", "", growth_rank[i]))
      gr_names[i]  = unique(as.character(csd_muts$cluster_id[csd_muts$log_rank == j]))
    }
    names(growth_rank) = gr_names
 
    csd_muts$table_rank = glue('{csd_muts$sel_lineage} : {csd_muts$sel_muts_all}')
    csd_select = csd_muts
    
    sel_mut_log_rank =unique(csd_muts$sel_mut_log_rank)
    sel_mut_lineages_log_rank = unique(csd_muts$sel_mut_lineages_log_rank)
    sel_lin_mut_log_rank = unique(csd_muts$sel_lineages_mut)
    log_rank = sub(".*,","", sel_mut_lineages_log_rank)
    sel_lineages = sub(",.*", "", sel_mut_lineages_log_rank) 
    sel_mutations = sub(",.*","", sel_lin_mut_log_rank)
    
    
    csd_lin = list()
    csd_names = NA
    
    for( j in 1:length(log_rank)) {
      csd_lin[[j]] = csd[csd$log_rank %in% log_rank[[j]],] 
      csd_lin[[j]]$sel_mutation = sel_mutations[[j]]
      csd_lin[[j]]$sel_mutation_pres = str_detect(csd_lin[[j]]$variants, sel_mutations[[j]])
        #ifelse(csd_lin[[j]]$sequence_name %in% csd_muts1[[sel_mutations[j]]]$sequence_name, 
         #                                     TRUE, FALSE)
      csd_names[j] = log_rank[[j]]
      
      
    }
    
    names(csd_lin) = csd_names
    
    csd_lin_all = dplyr::bind_rows(csd_lin)
    
    for(i in csd_lin_all$cluster_id){
    csd_lin_all$growth_rank[csd_lin_all$cluster_id == i] = growth_rank[names(growth_rank) == i]
    }
    csd_cluster = csd_lin_all
    perc_all = list() 
    
    for(j in 1:length(mutations)) { 
      
      perc_lin = list() 
      cluster_ids = unique(csd_muts1[[j]]$log_rank)
      
      for(i in 1:length(cluster_ids)) { 
        
        #percent of cluster that has mutation 
        perc_lin[[i]] = data.frame(table(csd_lin_all$sel_mutation_pres[csd_lin_all$sel_mutation == mutations[j] & 
                                                                         csd_lin_all$log_rank == cluster_ids[i]]))
        perc_lin[[i]]$percent = perc_lin[[i]]$Freq/sum(perc_lin[[i]]$Freq) * 100
        perc_lin[[i]]$log_rank = cluster_ids[i]
        perc_lin[[i]]$cluster_id = unique(csd_lin_all$cluster_id[csd_lin_all$log_rank == cluster_ids[i]])
        perc_lin[[i]]$sel_mutations = mutations[[j]]
        perc_lin[[i]] = perc_lin[[i]][perc_lin[[i]]$Var1 == TRUE,]
        
      }
      perc_lin = dplyr::bind_rows(perc_lin)
      perc_all[[j]] = perc_lin
      
    } 
    
    #names(perc_all) = mutations
    perc_all = bind_rows(perc_all)
    cluster_ids = log_rank
    #defining mutations 
    Ygr2 = Ygr1[Ygr1$log_rank %in% cluster_ids,]
    kp_nodes = as.numeric(unique(Ygr2$cluster_id))
    max_growth_rate = sprintf(as.numeric(max(Ygr2$logistic_growth_rate)), fmt = '%#.2f') 
    n_clusters_included = length(cluster_ids)
    cmts1 = cluster_muts(scanner_env = envs, nodes = kp_nodes, overlap_threshold = 90/100)
    cmts2 = cluster_muts_new(scanner_env = envs, nodes = kp_nodes, overlap_threshold = 90/100)
    cids2 = cluster_dels_all(scanner_env = envs, nodes = kp_nodes, overlap_threshold = 90/100)
    cids1 = cluster_dels(scanner_env = envs, nodes = kp_nodes, overlap_threshold = 90/100)
    names(cmts1) = sort(kp_nodes)
    names(cmts2) = sort(kp_nodes)
    names(cids1) = sort(kp_nodes)
    names(cids2) = sort(kp_nodes)
    #remove synsnps 
    e484k = list()
    s494p = list()
    q667h = list()
    e484k_all = list()
    s494p_all = list()
    l18f_all = list()
    q667h_all = list()
    for(i in 1:length(Ygr2$log_rank)) { 
      
      cmts1[[i]] = cmts1[[i]][!str_detect(cmts1[[i]],pattern="synSNP")]
      e484k[[i]] = cmts1[[i]][str_detect(cmts1[[i]],pattern="E484K")]
      q667h[[i]] = cmts1[[i]][str_detect(cmts1[[i]],pattern="Q677H")]
      s494p[[i]] = cmts1[[i]][str_detect(cmts1[[i]],pattern="S494P")]
      cmts3 = names(cmts2[[i]])
      e484k_all[[i]] = cmts3[str_detect(cmts3,pattern="E484K")]
      q667h_all[[i]] = cmts3[str_detect(cmts3,pattern="Q677H")]
      #l18f_all[[i]] = cmts3[str_detect(cmts3,pattern="Q677H")]
      s494p_all[[i]] = cmts3[str_detect(cmts3,pattern="S494P")]
    }
    names(e484k) = sort(kp_nodes)
    names(e484k_all) = sort(kp_nodes)
    names(q667h) = sort(kp_nodes)
    names(q667h_all) = sort(kp_nodes)
    names(s494p) = sort(kp_nodes)
    names(s494p_all) = sort(kp_nodes)
    
    dpa = list()
    dpa2 = list() 
    dpa3 = list()
    dpa4 = list() 
    for(i in 1:length(mutations)){
      
      for(j in 1:length(kp_nodes)){
        
        cmts3 = names(cmts2[[j]])
        dpa[[j]] = cmts3[str_detect(cmts3,pattern = mutations[[i]])] 
        dpa2[[j]] = cmts1[[j]][str_detect(cmts1[[j]],pattern = mutations[[i]])]
        
      }
      names(dpa) = sort(kp_nodes)
      names(dpa2) = sort(kp_nodes)
      dpa3[[i]] = dpa # present 
      dpa4[[i]] = dpa2 # defining
    }
    
    #setting up defining mutants, lineages, other lineages
    Ygr2 = list()
    
    for(j in 1:length(mutations)){
      
      Ygr2[[j]] = Ygr1[Ygr1$log_rank %in% csd_muts1[[j]]$log_rank,] 
      cluster_ids = unique(csd_muts1[[j]]$log_rank)
      
      for(i in cluster_ids) {
        
        lsr = as.integer(100*rev(sort(table(csd_lin_all$lineage[csd_lin_all$log_rank == i])/nrow(csd_lin_all[csd_lin_all$log_rank == i,]))))
        lsr = ifelse(lsr == 0, "<1", lsr)
        
        Ygr2[[j]]$percent[Ygr2[[j]]$log_rank == i] = perc_all$percent[perc_all$log_rank == i & perc_all$sel_mutations == mutations[j]]
        Ygr2[[j]]$lineages[Ygr2[[j]]$log_rank == i] = paste(unique(csd_lin_all$lineage_muts[csd_lin_all$log_rank == i]), collapse = ' ')
        Ygr2[[j]]$lineages_perc[Ygr2[[j]]$log_rank == i] = paste(glue('{names(rev(sort(table(csd_lin_all$lineage[csd_lin_all$log_rank == i]))))}({lsr})'), 
                                                             collapse = ' ')
        Ygr2[[j]]$table_rank[Ygr2[[j]]$log_rank == i] = unique(csd_muts$table_rank[csd_muts$log_rank == i])
        Ygr2[[j]]$sel_mutation[Ygr2[[j]]$log_rank == i] = mutations[j]
        Ygr2[[j]]$growth_rank[Ygr2[[j]]$log_rank == i] = unique(csd_muts$growth_rank[csd_muts$log_rank == i ])
        Ygr2[[j]]$table_rank[Ygr2[[j]]$log_rank == i] = unique(csd_muts$table_rank[csd_muts$log_rank == i ])
        Ygr2[[j]]$tips = NA
      }
    }
    
    
    for(i in sort(kp_nodes)) { 
      
      for(j in 1:length(mutations)) { 
        
        Ygr2[[j]]$defining_muts[Ygr2[[j]]$cluster_id == i] = paste(cmts1[[as.character(i)]], collapse = ' ') 
        Ygr2[[j]]$defining_dels[Ygr2[[j]]$cluster_id == i] = paste(cids1[[as.character(i)]], collapse = ' ')
        
        if(length(mutations) > 1) {
        Ygr2[[j]]$dpa[Ygr2[[j]]$cluster_id == i] = ifelse(!is_empty(dpa4[[j]][[as.character(i)]]) & 
                                                            !is_empty(dpa3[[j]][[as.character(i)]]), glue('{mutations[j]} Defining'), 
                                                          ifelse(is_empty(dpa4[[j]][[as.character(i)]]) &
                                                                   !is_empty(dpa3[[j]][[as.character(i)]]) ,glue('{mutations[j]} Present'), 
                                                                 glue('{mutations[j]} Absent')))
        
        } 
        
        if(length(mutations) == 1) {  
          
          Ygr2[[j]]$dpa[Ygr2[[j]]$cluster_id == i] = ifelse(!is_empty(dpa4[[j]][[as.character(i)]]) & 
                                                            !is_empty(dpa3[[j]][[as.character(i)]]), "Defining", 
                                                            ifelse(is_empty(dpa4[[j]][[as.character(i)]]) &
                                                                   !is_empty(dpa3[[j]][[as.character(i)]]) , "Present", 
                                                                   "Absent"))
          
      
          }
      }
      
    }
    
    if(!is.null(cut_off)){ 
      
      if(cut_off < 5){
        
        for(i in 1:length(mutations)) { 
          
          Ygr2[[i]] = Ygr2[[i]][Ygr2[[i]]$logistic_growth_rate >= cut_off,]
          
        } 
        
      } else { 
        
        for(i in 1:length(mutations)){ 
          
          Ygr2[[i]] = Ygr2[[i]][Ygr2[[i]]$log_rank <= cut_off,]

          
        }
        
        
      }
    }
    
   
    ygrl_all = list()
    
    
    for(i in 1:length(lineages)){
      
      NAs = rep(NA, length(mutations))
      ygrl = data.frame(group = NAs, table_rank = NAs, cluster_size_perc = NAs, lineages_perc = NAs, 
                        most_recent_tip = rep(ymd("2020-01-01"), length(mutations)), 
                        least_recent_tip = rep(ymd("2020-01-01"), length(mutations)), 
                        logistic_growth_rate = NAs, log_rank = NAs, 
                        dpa = NAs, defining_muts = NAs, defining_dels = NAs)
    
      
      for(j in 1:length(mutations)){ 
      
      perc_cluster = 100*nrow(csd_muts1[[j]][csd_muts1[[j]]$lineage == lineages[i],])/
        nrow(csd[csd$lineage == lineages[i],])
      if(length(mutations) > 1){ 
      ygrl$dpa[j] = ifelse(perc_cluster > 90, glue('{mutations[j]} Defining'), 
                        ifelse(perc_cluster > 0 & perc_cluster < 90,glue('{mutations[j]} Present'), 
                               glue('{mutations[j]} Absent')))
      } else { 
        ygrl$dpa[j] = ifelse(perc_cluster > 90, glue('Defining'), 
                           ifelse(perc_cluster > 0 & perc_cluster < 90,glue('Present'), 
                                  glue('Absent'))) 
        }
      
      perc_cluster = ifelse(perc_cluster >= 1, as.character(as.integer(perc_cluster)), 
                            ifelse(perc_cluster < 1 & perc_cluster > 0, "<1", 0))
      ygrl$group[j] = "Lineage"
      ygrl$table_rank[j] = lineages[i]
      ygrl$cluster_size_perc[j] = glue('{nrow(csd[csd$lineage == lineages[i],])}({perc_cluster})')
      ygrl$lineages_perc[j] = ""
      ygrl$most_recent_tip[[j]] = max(csd$most_recent_tip[csd$lineage == lineages[i]]) 
      ygrl$least_recent_tip[j] = min(csd$least_recent_tip[csd$lineage == lineages[i]])
      ygrl$logistic_growth_rate[j] = ""
      ygrl$log_rank[j] = ""
      ygrl$defining_muts = ""
      ygrl$defining_dels = ""
      
      }
      
      ygrl_all[[i]] = ygrl
      
    }
    
    ygrm_all = list()
    
    for(i in 1:length(mutations)){ 
    
      NAs = rep(NA, length(lineages))
      ygrm = data.frame(group = NAs, table_rank = NAs, cluster_size_perc = NAs, lineages_perc = NAs, 
                        most_recent_tip = rep(ymd("2020-01-01"), length(lineages)), 
                        least_recent_tip = rep(ymd("2020-01-01"), length(lineages)), 
                        logistic_growth_rate = NAs, log_rank = NAs, 
                        dpa = NAs, defining_muts = NAs, defining_dels = NAs)
    
    for(j in 1:length(lineages)){ 
    
      perc_cluster = 100 * nrow(csd_muts1[[i]][csd_muts1[[i]]$lineage == lineages[j],])/
        nrow( csd %>% filter(str_detect(variants, mutations[[i]])))
      
      if(length(lineages) > 1){ 
        ygrm$dpa[j] = ifelse(perc_cluster >= 90, glue('{lineages[j]} Defining'), 
                           ifelse(perc_cluster > 0 & perc_cluster < 90, glue('{lineages[j]} Present'), 
                                  glue('{lineages[j]} Absent')))
        
      } else { 
        
       ygrm$dpa[j] = ifelse(perc_cluster >= 90, glue('Defining'), 
                             ifelse(perc_cluster > 0 & perc_cluster < 90, glue('Present'), 
                                    glue('Absent')))
       
      }
      perc_cluster = ifelse(perc_cluster >= 1, as.character(as.integer(perc_cluster)), 
                            ifelse(perc_cluster < 1 & perc_cluster > 0, "<1", 0))
      ygrm$cluster_size_perc[[j]] = glue('{nrow( csd %>% filter(str_detect(variants, mutations[[i]])))}({perc_cluster})')
      ygrm$most_recent_tip[[j]] = max(csd$most_recent_tip[csd$lineage == lineages[j]]) 
      ygrm$least_recent_tip[[j]] = min(csd$least_recent_tip[csd$lineage == lineages[j]])
      
    }  
    
    ygrm$group = "Mutation"
    ygrm$table_rank = mutations[i]
    ygrm$lineages_perc = ""
    ygrm$logistic_growth_rate = ""
    ygrm$log_rank = ""
    ygrm$defining_muts = ""
    ygrm$defining_dels = ""
    
    ygrm_all[[i]] = ygrm
    
    }
    
    ygrl_all = dplyr::bind_rows(ygrl_all, ygrm_all)
    
    Ygr2 = dplyr::bind_rows(Ygr2)
    Ygr2$logistic_growth_rate = sprintf(as.numeric(Ygr2$logistic_growth_rate), fmt = '%#.2f')   
    Ygr2$percent = as.integer(Ygr2$percent) 
    #Ygr2 = Ygr2[Ygr2$percent > prop_cluster,]
    Ygr2$cluster_size_perc = glue('{Ygr2$cluster_size}({Ygr2$percent})')
    Ygr2$group = "Selection"
    #Ygr2$lineages_perc = Ygr2$lineages2 #glue('{Ygr2$lineages}({Ygr2$percent})')
    #Ygr2$lineages2 = NULL 
    Ygr2$log_rank = as.character(Ygr2$log_rank)
    
    Ygr2 = dplyr::bind_rows(ygrl_all, Ygr2)
    Ygr2_output = Ygr2 %>% select(group, table_rank, cluster_size_perc, lineages_perc,  most_recent_tip, least_recent_tip, logistic_growth_rate, log_rank, dpa,  defining_dels, defining_muts)
    Ygr2_output = Ygr2_output[order(Ygr2_output$group),]
    Ygr2_output = Ygr2_output[,-1]
    colnames(Ygr2_output) = c("Selection", "Size (%)", "Lineage (%)", "Most recent tip", 
                              "Least recent tip", "Log growth rate", "Growth rate rank", 
                              "Mutation status", "Defining deletions", "Defining mutations")
    rownames(Ygr2_output) <- c()
    
    
    
  }


saveRDS(growth_rank, glue('{path_to_scanner}growth_rank.rds'))
saveRDS(Ygr2, glue('{path_to_scanner}Ygr2_{max_date}_{min_date}.RDS'))
kp_nodes = kp_nodes[kp_nodes %in% Ygr2$cluster_id]
saveRDS(kp_nodes, glue('{path_to_scanner}kp_nodes.RDS'))

Ygr2_out = list(Ygr2_output, csd_select, csd_cluster, kp_nodes, growth_rank, cmts1, cmts2, Ygr2)

}





