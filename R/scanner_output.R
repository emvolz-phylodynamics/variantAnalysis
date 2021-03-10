#' Output a markdown for scanner results, including full written report on selected inclusion criteria of interest 
#' 
#' Currently includes cumulative sum of sampled sequences, logistic growth rate and frequency in log odds over time,
#' coordinates of logistic growth rates, lineages/clusters of interest, samples dates, 
#' growth rates over time for lineages of interest, proportion of sequences with mutations of interest. 
#' 
#' Age proportions in cluster can be added (if scanner outputs stored locally and linked with PHE, to be added in...)
#' 
#' TO ADD: treestructure figure linking growth rates to node number and p-value on tree for each cluster of interest
#' Keep random sample of size N from tree with sequences from clusters of interest, node values on tree
#' Clusters coloured using Dark2 
#' 
#' @param lineages Lineage of interest, if NULL not used for inclusion criteria
#' @param mutations Mutations of interest, if NULL not used for inclusion criteria
#' @param cut_off Statistical cut off of interest, if NULL not used for inclusion criteria. If growth rate below 
#' 5 selected, will use this as a Growth rates > x cut off. If cut_off >= 5 selected, will return the clusters 
#' with top x growth rates. 
#' @param path_to_scanner Path to scanner run outputs and scanner metadata (scanner.rds, scanner.csv & scanner_env)
#' @param path_to_data Path to COG tree and COG tree metadata (sample_date, )
#' @param path_to_metadata Path to sequence metadata (age, gender) e.g. PHE data, if NULL no age plots included
#' @param log_growth_rate_cut_off Logistic growth rate cut-off we use as exclusion criteria 
#' @param defining_mutations_cut_off Minimm proportion (%) of sequences in clusters of interest and maximum proportion (%) of 
#' of sequences in comparison general population cluster that must have a specific mutation for it to be considered
#' a defining mutation of a cluster of interest
#' @param prop_mutation Proportion of sequences (%) with mutations of interest across selected clusters
#' @param max_date Most recent tip sample date from scanner run 
#' @param min_date Least recent tip sample date from scanner run 
#' @param path_to_save Path to directory to save markdown in, if NULL will save in current working directory  
#' @export 

scanner_output = function(lineages = c("A.23.1", "B.1.525", "B.1.351"), mutations = NULL, cut_off = NULL, 
         path_to_scanner = "C:/Users/oboyd/Documents/COGUK_server/min_10/3_month/", 
         scanner_run_date = "2020-03-06",
         path_to_data = "C:/Users/oboyd/Documents/COGUK_server/", path_to_metadata = NULL, 
         log_growth_rate_cut_off = 0.0, defining_mutations_cut_off = 90, 
         prop_mutations = 25, path_to_save = NULL, min_size = 10, max_size = 20000,
         max_date = "2021-02-19", min_date = "2020-11-01", 
         include_pillar1 = FALSE, scanner_run_Date = "2021-02-16", 
         generation_time = 6.5) { 
  
  library(glue)
  
  if(!is.null(path_to_save)) { 
    
  rmarkdown::render('master_scanner.Rmd' 
                    ,output_file = glue('scanner-{max_date}-{min_date}-{ifelse(!is.null(lineages),paste0(paste(lineages, collapse = "-"),"-"),"")}{ifelse(!is.null(mutations),paste0(paste(mutations, collapse = "-"),"-"),"")}{ifelse(!is.null(cut_off),paste0(paste(cut_off, collapse = "-"),"-"),"")}report') 
                    ,output_dir = path_to_save
                    ) 
  } else { 
    
    rmarkdown::render('master_scanner.Rmd'
                      ,output_file = glue('scanner-{max_date}-{min_date}-{ifelse(!is.null(lineages),paste0(paste(lineages, collapse = "-"),"-"),"")}{ifelse(!is.null(mutations),paste0(paste(mutations, collapse = "-"),"-"),"")}{ifelse(!is.null(cut_off),paste0(paste(cut_off, collapse = "-"),"-"),"")}report')
                      )
  }
  
}

#' to run: scanner_output()

