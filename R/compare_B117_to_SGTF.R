# 
#  _   _                    ____   ____ _____ _____ 
# | \ | | ___  __   _____  / ___| / ___|_   _|  ___|
# |  \| |/ _ \ \ \ / / __| \___ \| |  _  | | | |_   
# | |\  |  __/  \ V /\__ \  ___) | |_| | | | |  _|  
# |_| \_|\___|   \_/ |___/ |____/ \____| |_| |_|    
# 

require(hrbrthemes)
require(scales)

# this is copied and pasted from compare_lineages_d2.R!
B.1.1.7_matched_notB.1.1.7_deduped_meanrate0.0005 = unpack_two_Lineages(ofn = "C:/Users/lilyl/OneDrive/Documents/N501Y_analyses/results/2021-01-18/d1-tN-sampler1_B.1.1.7_vs_matchSample_notB.1.1.7_meanrate5e-04-dedup.rds",
                                                                        Lineage_main = 'B.1.1.7',
                                                                        Lineage_matched = 'Control',
                                                                        dedup = "-dedup",
                                                                        meanrate = 0.0005)




mleskystp <- readRDS("~/N501Y_analyses/dropouts/mleskystp.rds")

sgss_stp_new_43_52_weeks <- readRDS("~/N501Y_analyses/dropouts/sgss_stp_new_43_52_weeks.rds")

# plot(mleskystp[[1]][[2]]$ne)

sgss_stp_new_43_52_weeks$true_negative = sgss_stp_new_43_52_weeks$sgss_s_negative * sgss_stp_new_43_52_weeks$total_cases /
  (sgss_stp_new_43_52_weeks$sgss_s_negative + sgss_stp_new_43_52_weeks$sgss_s_positive)


sgss_stp_new_43_52_weeks$true_negative_corrected = sgss_stp_new_43_52_weeks$sgss_s_negative_corrected * sgss_stp_new_43_52_weeks$total_cases /
  (sgss_stp_new_43_52_weeks$sgss_s_negative_corrected + sgss_stp_new_43_52_weeks$sgss_s_positive_corrected)





pldf <- as.data.frame(do.call(rbind, lapply(unique(sgss_stp_new_43_52_weeks$epiweek), function(week) {
  if(nrow(sgss_stp_new_43_52_weeks[sgss_stp_new_43_52_weeks$epiweek == week, ]) == 42) {# checking there are no duplicates; there are 42 STPs
    total_S_neg = sum(sgss_stp_new_43_52_weeks[sgss_stp_new_43_52_weeks$epiweek == week, "sgss_s_negative_corrected"])
    
    ne = median( B.1.1.7_matched_notB.1.1.7_deduped_meanrate0.0005$q_ne[which(lubridate::epiweek(lubridate::date_decimal(B.1.1.7_matched_notB.1.1.7_deduped_meanrate0.0005$time)) == week), ][,"y"])
    neub = median( B.1.1.7_matched_notB.1.1.7_deduped_meanrate0.0005$q_ne[which(lubridate::epiweek(lubridate::date_decimal(B.1.1.7_matched_notB.1.1.7_deduped_meanrate0.0005$time)) == week), ][,"yub"])
    nelb = median( B.1.1.7_matched_notB.1.1.7_deduped_meanrate0.0005$q_ne[which(lubridate::epiweek(lubridate::date_decimal(B.1.1.7_matched_notB.1.1.7_deduped_meanrate0.0005$time)) == week), ][,"ylb"])
    
    
    return(c(week = week, total_S_neg = total_S_neg, ne = ne, neub = neub, nelb = nelb))
  }
}
)
)
)


pl = ggplot(pldf) + geom_point(aes(x = ne, y = total_S_neg), shape = 15) + 
  geom_errorbarh(aes(xmin = nelb, xmax = neub, y = total_S_neg)) +
  # scale_x_log10() + scale_y_log10() +
  scale_y_continuous(
    trans = "log10",
    breaks = function(x) {
      brks <- extended_breaks(Q = c(1, 5))(log10(x))
      10^(brks[brks %% 1 == 0])
    },
    labels = math_format(format = log10)
  ) +
  scale_x_continuous(
    trans = "log10",
    breaks = function(x) {
      brks <- extended_breaks(Q = c(1, 5))(log10(x))
      10^(brks[brks %% 1 == 0])
    },
    labels = math_format(format = log10)
  )+
  theme_bw() + labs(x = "Effective population size B.1.1.7", y = "TPR-adjusted SGTF case numbers") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=14,face="bold"), panel.grid.minor = element_blank())+ annotation_logticks()
pl

ggsave( plot = pl, file = "TPR-adjusted_SGTF_vs_Ne_mlesky.pdf", width = 8, height = 8 )



### just to see what it looks like with raw S-

pldf <- as.data.frame(do.call(rbind, lapply(unique(sgss_stp_new_43_52_weeks$epiweek), function(week) {
  if(nrow(sgss_stp_new_43_52_weeks[sgss_stp_new_43_52_weeks$epiweek == week, ]) == 42) {# checking there are no duplicates; there are 42 STPs
    total_S_neg = sum(sgss_stp_new_43_52_weeks[sgss_stp_new_43_52_weeks$epiweek == week, "sgss_s_negative_corrected"])
    total_S_neg_raw = sum(sgss_stp_new_43_52_weeks[sgss_stp_new_43_52_weeks$epiweek == week, "sgss_s_negative"])
    
    ne = median( B.1.1.7_matched_notB.1.1.7_deduped_meanrate0.0005$q_ne[which(lubridate::epiweek(lubridate::date_decimal(B.1.1.7_matched_notB.1.1.7_deduped_meanrate0.0005$time)) == week), ][,"y"])
    neub = median( B.1.1.7_matched_notB.1.1.7_deduped_meanrate0.0005$q_ne[which(lubridate::epiweek(lubridate::date_decimal(B.1.1.7_matched_notB.1.1.7_deduped_meanrate0.0005$time)) == week), ][,"yub"])
    nelb = median( B.1.1.7_matched_notB.1.1.7_deduped_meanrate0.0005$q_ne[which(lubridate::epiweek(lubridate::date_decimal(B.1.1.7_matched_notB.1.1.7_deduped_meanrate0.0005$time)) == week), ][,"ylb"])
    
    
    return(c(week = week, TPR_adjusted = total_S_neg, Raw = total_S_neg_raw, ne = ne, neub = neub, nelb = nelb))
  }
}
)
)
)
names(pldf) = c("Week", "TPR-adjusted", "Raw", "ne", "neub" ,     "nelb")
pldf = reshape2::melt(pldf, id.vars = names(pldf)[!names(pldf) %in% c("Raw" , "TPR-adjusted"  )])

tmp_unadj = pldf[pldf$variable == "Raw", ]
tmp_adj = pldf[pldf$variable == "TPR-adjusted", ]

my.formula_unadj <- log(tmp_unadj$value) ~ log(tmp_unadj$ne)
my.formula_adj <- log(tmp_adj$value) ~ log(tmp_adj$ne)




pl_2 = ggplot(pldf,aes(x = ne, y = value, col = variable, fill = variable)) + geom_point( shape = 15, size = 2)  + 
  geom_errorbarh(aes(xmin = nelb, xmax = neub, y = value, col = variable)) + 
  # scale_x_log10() + 
  scale_y_continuous(
    trans = "log10",
    breaks = function(x) {
      brks <- extended_breaks(Q = c(1, 5))(log10(x))
      10^(brks[brks %% 1 == 0])
    },
    labels = math_format(format = log10)
  ) + 
  scale_x_continuous(
    trans = "log10",
    breaks = function(x) {
      brks <- extended_breaks(Q = c(1, 5))(log10(x))
      10^(brks[brks %% 1 == 0])
    },
    labels = math_format(format = log10)
  )+    
  stat_smooth(method = "lm")+
  
  theme_bw() + labs(x = "Effective population size B.1.1.7", y = "SGTF case numbers", col = "SGTF", fill = "SGTF") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=14,face="bold"), legend.position = c(0.7,0.3), legend.title = element_text(size=14), 
        legend.text =element_text(size=12) ,legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'),
        plot.margin = margin(10, 10, 10, 10), panel.grid.minor = element_blank())+ 
  annotation_logticks() +
  ggpmisc::stat_poly_eq(data = tmp_unadj, formula = my.formula_unadj, size = 5,
                        mapping = aes(log(ne), log(value),  col = variable, fill = variable,
                                      label = paste( ..adj.rr.label.., sep = "~~~")),  parse = TRUE)+  
  ggpmisc::stat_poly_eq(data = tmp_adj, label.y = 0.85, size = 5,
                        formula = my.formula_adj, mapping = aes(log(ne), log(value), col = variable, fill = variable,
                                                                label = paste( ..adj.rr.label.., sep = "~~~")),  parse = TRUE)

pl_2


ggsave( plot = pl_2, file = "TPR-adjusted_SGTF_vs_Ne_mlesky_raw_vs_adjusted.pdf", width = 8, height = 8 )


