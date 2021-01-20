library(sarscov2) 
library(ape) 
library(lubridate)
library(treedater) 
require(ggplot2)
require(grid)
require(gridExtra)
require(ggtree)
library(alakazam)
require(stringi)

source("~/R/erik_skygrowth_functions.R")
source("~/R/functions.R")
source("~/R/sampling.R")
source("~/R/dating_and_mlesky_functions.R")


root_dir = '/cephfs/covid/bham/climb-covid19-volze'
regenerate_ML_tree = F
lineage = 'B.1.1.7' 
weightsfn = '/cephfs/covid/bham/climb-covid19-volze/b0-weightsdf-2021-01-04.csv'
mindate = as.Date( '2020-10-15' ) 
maxdate = as.Date( Sys.Date() - 12 )
ns = c( 1000 )
prop_stratified = .25
deduplicate = TRUE 




lsoa_to_stp = read.csv(paste0(root_dir, '/lsoa_to_stp.csv'))
pcdf = read.csv( '/cephfs/covid/bham/climb-covid19-geidelbergl/md/postcode_to_la.csv', header=TRUE, stringsAs=FALSE)
md = read.csv(paste0(root_dir, '/phylolatest/alignments/', 'cog_2021-01-09_metadata.csv'))


mintime <- decimal_date( mindate ) 
maxtime <- decimal_date( maxdate )

wdf = na.omit( read.csv( weightsfn , stringsAs=FALSE )  ) 


civetfn =  list.files(  paste0(root_dir, '/phylolatest/civet/' ), patt = 'cog_global_[0-9\\-]+_metadata.csv', full.names=TRUE) #'../phylolatest/civet/cog_global_2020-12-01_metadata.csv'
civmd = read.csv( civetfn , stringsAs=FALSE , header=TRUE )
civmd$central_sample_id <-  sapply( strsplit( civmd$sequence_name , split='/' ) , '[', 2 ) # for linkage 
civmd$sample_date <- as.Date( civmd$sample_date )
civmd$sample_time <- decimal_date( civmd$sample_date ) 

# filter dates 
civmd <- civmd[ civmd$sample_time >= mintime , ]
civmd <- civmd[ civmd$sample_time <= maxtime , ]

# filter lineage 
civmd <- civmd[ civmd$lineage == lineage , ]

# combine
s <- merge( civmd, wdf, by = 'central_sample_id' )

# exclude p1 
lhls <- c( 'MILK', 'ALDP', 'QEUH', 'CAMC')
s$lighthouse = FALSE
for (lh in lhls ){
  s$lighthouse[ grepl(s$central_sample_id, patt = lh) ] <- TRUE
}
s <- s [ s$lighthouse , ]
s$epiweek <- epiweek(  s$sample_date )


stp_list = lapply(names(table(lsoa_to_stp$STP19NM)), function(region) {
  postcodes <- pcdf[ pcdf$LTLA19NM %in%   lsoa_to_stp[lsoa_to_stp$STP19NM == region, "LAD19NM"], "postcode"]
  seq_in_region = sapply(strsplit(as.character(md[ md$outer_postcode %in% postcodes , "sequence_name"]), '/'), '[[', 2)
  
  
  seq_in_region = s[s$central_sample_id %in% seq_in_region, 'central_sample_id']
  
  if(length(seq_in_region) >= 100) {
    return(list(STP = region,
                sequences = seq_in_region))
  } else {
    return(NULL)
  }
}

)

stp_list = stp_list[!unlist(lapply(stp_list, is.null))]


# load tree and deduplicate 
tr0 = read.tree( list.files( paste0(root_dir, '/phylolatest/trees'), patt = '.*newick', full=T ))
tr0$tip.label <-  sapply( strsplit( tr0$tip.label, split='/' ), '[', 2 )

lapply(stp_list, function(x) {
  if(!is.null(x)) {
    
    tres <- list(keep.tip(tr0, x$sequences),
                 keep.tip(tr0, x$sequences))
    class( tres ) <- 'multiPhylo' 
    x$STP = gsub(' ','-',x$STP)
    write.tree( tres, file = glue('sampler1_{lineage}_{Sys.Date()}_n={length(x$sequences)}{x$STP}.nwk')  )}
  
  
  
})

# tr1 <- keep.tip( tr0, s$central_sample_id )


bootrep_2 <- function(mltr, ...)
{
  td = datetree( mltr , ... ) 
  res = diff( range( epiweek( date_decimal( td$sts )  )  )  ) + 1
  res <- res * 2
  sg = mlskygrid( td, tau = NULL, tau_lower=.001, tau_upper = 10 , sampleTimes = td$sts , res = 12, ncpu = 6)
  sg
}

out = lapply(stp_list, function(x) {
  if(!is.null(x)) {
    tres <- list(keep.tip(tr0, x$sequences))
    class( tres ) <- 'multiPhylo' 
    meanrate = .001
    taxis <- decimal_date( seq( as.Date( mindate) , as.Date(maxdate)  , by = 1) )
    return(list(loc = x$STP, bootrep_2(tres[[1]])))
  }
  
  
})

saveRDS(out, 'mleskystp.rds')

# 
#   ____ _____ _____   ____    _  _____  _    _____ ____      _    __  __ _____   _____ ____   ___  __  __    ____ _     ___ __  __ ____  
#  / ___| ____|_   _| |  _ \  / \|_   _|/ \  |  ___|  _ \    / \  |  \/  | ____| |  ___|  _ \ / _ \|  \/  |  / ___| |   |_ _|  \/  | __ ) 
# | |  _|  _|   | |   | | | |/ _ \ | | / _ \ | |_  | |_) |  / _ \ | |\/| |  _|   | |_  | |_) | | | | |\/| | | |   | |    | || |\/| |  _ \ 
# | |_| | |___  | |   | |_| / ___ \| |/ ___ \|  _| |  _ <  / ___ \| |  | | |___  |  _| |  _ <| |_| | |  | | | |___| |___ | || |  | | |_) |
#  \____|_____| |_|   |____/_/   \_\_/_/   \_\_|   |_| \_\/_/   \_\_|  |_|_____| |_|   |_| \_\\___/|_|  |_|  \____|_____|___|_|  |_|____/ 
#   
#   > 
##############

require(hrbrthemes)

mleskystp <- readRDS("~/N501Y_analyses/dropouts/mleskystp.rds")

sgss_stp_new_43_52_weeks <- readRDS("~/N501Y_analyses/dropouts/sgss_stp_new_43_52_weeks.rds")

# plot(mleskystp[[1]][[2]]$ne)

sgss_stp_new_43_52_weeks$true_negative = sgss_stp_new_43_52_weeks$sgss_s_negative * sgss_stp_new_43_52_weeks$total_cases /
  (sgss_stp_new_43_52_weeks$sgss_s_negative + sgss_stp_new_43_52_weeks$sgss_s_positive)


sgss_stp_new_43_52_weeks$true_negative_corrected = sgss_stp_new_43_52_weeks$sgss_s_negative_corrected * sgss_stp_new_43_52_weeks$total_cases /
  (sgss_stp_new_43_52_weeks$sgss_s_negative_corrected + sgss_stp_new_43_52_weeks$sgss_s_positive_corrected)



# sgss_stp_new_43_52_weeks$true_negative = sgss_stp_new_43_52_weeks$true_negative_corrected



res = data.frame(do.call(rbind, lapply(mleskystp, function(x) {
  cbind(area = x[[1]], epiweek= lubridate::epiweek(lubridate::date_decimal(x[[2]]$time)) , x[[2]]$ne_ci)
})))

res = merge(res, sgss_stp_new_43_52_weeks)

res$epiweek = as.numeric(as.character(res$epiweek))
res$nelb = as.numeric(as.character(res$nelb))
res$ne = as.numeric(as.character(res$ne))
res$neub = as.numeric(as.character(res$neub))

slope = 643.04
intercept = 4306.91


# trajectories


ggplot(res, aes(x = epiweek, y = ne)) + geom_line() + facet_wrap(~area, scales = "free") + theme_bw() +
  geom_line(aes(x = epiweek , y = neub), linetype = 'dashed')+
  geom_line(aes(x = epiweek , y = nelb), linetype = 'dashed') + 
  geom_point(aes(x = epiweek, y = true_negative_corrected/slope ), col = 'red') + scale_y_continuous(
    
    # Features of the first axis
    name = "Ne",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*slope, name="True negative corrected SGTF")
  ) +  theme_ipsum() +
  
  theme(
    axis.title.y = element_text(color = 'black', size=33),
    axis.title.y.right = element_text(color = 'red', size=33)
  ) +
  
  ggtitle("Ne vs SGTF")


# res = data.frame(do.call(rbind, lapply(mleskystp, function(x) {
#   c(x[[1]], tail(x[[2]]$ne, 1))
# })),
# sgss_stp_new_43_52_weeks[sgss_stp_new_43_52_weeks$area %in% unlist(lapply(mleskystp, '[[', 1)) & sgss_stp_new_43_52_weeks$epiweek == 52, "sgss_s_negative"]
# )
# 
# names(res) = c("STP", "median_ne_week_52", "sgss_s_negative_week_52")
# res$median_ne_week_52 = as.numeric(as.character(res$median_ne_week_52))
# res$sgss_s_negative_week_52 = as.numeric(res$sgss_s_negative_week_52)
# plot(res$median_ne_week_52, res$sgss_s_negative_week_52)
# 
# 
# fit1 <- lm(sgss_s_negative_week_52 ~ median_ne_week_52, data = res)
# summary(fit1)
# 
# require(ggplot2)
# library(ggrepel)
# 
# ggplot(res, aes(median_ne_week_52, sgss_s_negative_week_52, label = STP)) +  theme_bw()+ geom_point(size = 3)+
#   stat_smooth(method = "lm", col = "red")+
#   geom_label_repel(alpha = 0.2)
# 
# 
# 
# 
# 
# 
# res = data.frame(do.call(rbind, lapply(mleskystp, function(x) {
#   c(x[[1]], tail(x[[2]]$ne, 1))
# })),
# sgss_stp_new_43_52_weeks[sgss_stp_new_43_52_weeks$area %in% unlist(lapply(mleskystp, '[[', 1)) & sgss_stp_new_43_52_weeks$epiweek == 52, "sgss_s_negative_corrected"]
# )
# 
# names(res) = c("STP", "median_ne_week_52", "sgss_s_negative_corrected_week_52")
# res$median_ne_week_52 = as.numeric(as.character(res$median_ne_week_52))
# res$sgss_s_negative_corrected_week_52 = as.numeric(res$sgss_s_negative_corrected_week_52)
# plot(res$median_ne_week_52, res$sgss_s_negative_corrected_week_52)
# 
# 
# fit1 <- lm(sgss_s_negative_corrected_week_52 ~ median_ne_week_52, data = res)
# summary(fit1)
# 
# require(ggplot2)
# library(ggrepel)
# 
# ggplot(res, aes(median_ne_week_52, sgss_s_negative_corrected_week_52, label = STP)) +  theme_bw()+ geom_point(size = 3)+
#   stat_smooth(method = "lm", col = "red")+
#   geom_label_repel(alpha = 0.2)


the_week = 43
# true_negative
# x = mleskystp[[6]]
res = data.frame(do.call(rbind, lapply(mleskystp, function(x) {
  c(x[[1]],  x[[2]]$ne_ci[which(lubridate::epiweek(lubridate::date_decimal(x[[2]]$time)) == the_week)[1],] )
})),
sgss_stp_new_43_52_weeks[sgss_stp_new_43_52_weeks$area %in% unlist(lapply(mleskystp, '[[', 1)) & sgss_stp_new_43_52_weeks$epiweek == the_week, "true_negative_corrected"]
)

names(res) = c("STP", "lowerCI_ne_week_52","median_ne_week_52","upperCI_ne_week_52", "true_negative_corrected")
res$lowerCI_ne_week_52 = as.numeric(as.character(res$lowerCI_ne_week_52))
res$upperCI_ne_week_52 = as.numeric(as.character(res$upperCI_ne_week_52))
res$median_ne_week_52 = as.numeric(as.character(res$median_ne_week_52))
res$true_negative_corrected = as.numeric(res$true_negative_corrected)


fit1 <- lm(true_negative_corrected ~ median_ne_week_52, data = res)
summary(fit1)

require(ggplot2)
library(ggrepel)
require(hrbrthemes)
# require(EpiEstim)

ggplot(res, aes(median_ne_week_52, true_negative_corrected, label = STP)) +  theme_bw()+ geom_point(size = 3)+
  stat_smooth(method = "lm", col = "red")+
  geom_label_repel(alpha = 0.2) + geom_errorbarh(aes(xmax = upperCI_ne_week_52, xmin = lowerCI_ne_week_52))







###


res = data.frame(do.call(rbind, lapply(mleskystp, function(x) {
  cbind(area = x[[1]], epiweek= lubridate::epiweek(lubridate::date_decimal(x[[2]]$time)) , x[[2]]$ne_ci)
})))

res = merge(res, sgss_stp_new_43_52_weeks)

res$epiweek = as.numeric(as.character(res$epiweek))
res$nelb = as.numeric(as.character(res$nelb))
res$ne = as.numeric(as.character(res$ne))
res$neub = as.numeric(as.character(res$neub))


# lm for each week
#   week = unique(res$epiweek)[[1]]
r2 = lapply(unique(res$epiweek), function(week) {
  
  tmp = res[res$epiweek == week, ]
  tmp = tmp[!duplicated(tmp$area) , ]
  fit <- lm(log(true_negative) ~ log(ne), data = tmp)
  
  pl <- ggplot(tmp, aes(ne, true_negative_corrected, label = area )) +  theme_bw()+ geom_point(size = 3)+
    stat_smooth(method = "lm", col = "blue", fill = 'blue')+
    geom_label_repel(alpha = 0.2) + geom_errorbarh(aes(xmax = neub, xmin = nelb)) + 
    facet_wrap(~epiweek) + 
    theme(strip.text.x = element_text(size = 20)) + 
    geom_point(size = 3, aes(ne, true_negative ), col = 'red')+
    stat_smooth(aes(ne, true_negative ), method = "lm", col = "red", fill = 'red')
  
  fit2 <- lm(log(true_negative_corrected) ~ log(ne), data = tmp)
  
  
  return(list(pl = pl, r2_raw = summary(fit)$adj.r.squared, r2_corrected = summary(fit2)$adj.r.squared, 
              p_raw = summary(fit)$coefficients[,4][2], p_corrected = summary(fit2)$coefficients[,4][2]  ))
}
)




plist = lapply(r2, '[[', 1)
library(gridExtra)
n <- length(plist)
nCol <- floor(sqrt(n))
do.call("grid.arrange", c(plist, ncol=nCol))










ggplot(data.frame(epiweek = unique(res$epiweek), r2_raw = unlist(lapply(r2, '[[', 2)), r2_corrected = unlist(lapply(r2, '[[', 3)))) +
  geom_point(aes(x = epiweek, y = r2_raw), col = 'red', size = 3)+
  geom_point(aes(x = epiweek, y = r2_corrected), col = 'blue', size = 3) + theme_bw() + labs(y = 'adjusted r^2 of linear regression log(ne) vs log(SGTF)') +
  annotate('text', x = 45, y = 0.7, label = 'raw SGSS S-', size = 9, col = 'red')+
  annotate('text', x = 45, y = 0.67, label = 'corrected SGSS S-', size = 9, col = 'blue') + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  scale_x_continuous(breaks = c(unique(res$epiweek)))



ggplot(data.frame(epiweek = unique(res$epiweek), p_raw = unlist(lapply(r2, '[[', 4)), p_corrected = unlist(lapply(r2, '[[', 5)))) +
  geom_point(aes(x = epiweek, y = log(p_raw)), col = 'red', size = 3)+
  geom_point(aes(x = epiweek, y = log(p_corrected)), col = 'blue', size = 3) + theme_bw() + labs(y = 'p-value (log scale) of lm of linear regresion log(ne) vs log(SGTF)')+
  annotate('text', x = 45, y = log(0.002), label = 'raw SGSS S-', size = 9, col = 'red')+
  annotate('text', x = 45, y = log(0.001), label = 'corrected SGSS S-', size = 9, col = 'blue')+ 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  scale_x_continuous(breaks = c(unique(res$epiweek)))


res2 = merge(res, data.frame(
  area = unlist(lapply(mleskystp, function(x) x[[1]])),
  nseq = unlist(lapply(mleskystp, function(x) length(x[[2]]$tre$sts))))
)
res3 = res2[res2$epiweek == 45, ]

fit3 = lm(ne ~ nseq , res3)
summary(fit3)

ggplot(res2, aes(x = nseq, y = ne )) + geom_point() + theme_bw() + facet_wrap(~epiweek, scales = 'free') +
  labs(x = 'Number of sequences', y = 'Estimate Ne')+ 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))
# WEEK 46??





# LOG-LOG!

# lm for each week
#   week = unique(res$epiweek)[[1]]
r2 = lapply(unique(res$epiweek), function(week) {
  
  tmp = res[res$epiweek == week, ]
  tmp = tmp[!duplicated(tmp$area) , ]
  # fit <- lm(log(true_negative_corrected) ~ log(ne), data = tmp)
  
  tmp$epiweek = paste('Week', tmp$epiweek)
  
  my.formula <- log(tmp$true_negative_corrected) ~ log(tmp$ne)
  
  
  pl <- ggplot(tmp, aes(log(ne), log(true_negative_corrected), label = area )) +  theme_bw()+ geom_point(size = 3)+
    stat_smooth(method = "lm", col = "red")+
    # geom_label_repel(alpha = 0.2) + geom_errorbarh(aes(xmax = log(neub), xmin = log(nelb))) + 
    # geom_label_repel(data = tmp[tmp$area == "Kent and Medway",], alpha = 0.7, arrow = NULL, nudge_x = log(0.496),
    #                  segment.color = NA) + 
    geom_errorbarh(aes(xmax = log(neub), xmin = log(nelb))) + 
    facet_wrap(~epiweek) + 
    theme(strip.text.x = element_text(size = 20)) + 
    labs(y = "TPR-adjusted SGTF (log scale)", x = "Effective population size (log scale)")+ 
    theme(axis.text=element_text(size=16),
          axis.title=element_text(size=14,face="bold"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()
    ) +
    ggpmisc::stat_poly_eq(data = tmp, formula = my.formula,
                          mapping = aes(log(ne), log(true_negative_corrected), label = paste( ..adj.rr.label.., sep = "~~~")),  parse = TRUE)
  
  
  fit2 <- lm(log(true_negative_corrected) ~ log(ne), data = tmp)
  
  
  return(list(pl = pl,# r2_raw = summary(fit)$adj.r.squared, 
              r2_corrected = summary(fit2)$adj.r.squared, 
              # p_raw = summary(fit)$coefficients[,4][2], 
              p_corrected = summary(fit2)$coefficients[,4][2]  ))
}
)




plist = lapply(r2, '[[', 1)
# library(gridExtra)
# n <- length(plist)
# nCol <- floor(sqrt(n))
# do.call("grid.arrange", c(plist, ncol=nCol))


library(grid)
library(gridExtra)

pl = cowplot::plot_grid(plotlist = plist, align = "v")


y.grob <- textGrob("TPR-adjusted SGTF (log scale)", 
                   gp=gpar(fontface="bold", fontsize=20), rot=90)

x.grob <- textGrob("Effective population size (log scale)", 
                   gp=gpar(fontface="bold",  fontsize=20))


grid.arrange(arrangeGrob(pl, left = y.grob, bottom = x.grob))


# save plot plot

pdf(file = 'results/SGTF_by_Ne_SG.pdf', width = 12, height = 12)
grid.arrange(arrangeGrob(pl, left = y.grob, bottom = x.grob))
dev.off()


png(file = 'results/SGTF_by_Ne_SG.png', width = 1000, height = 1000)
grid.arrange(arrangeGrob(pl, left = y.grob, bottom = x.grob))
dev.off()





#                                      _                         _  _           _           _                                  
#   ___ ___  _ __ ___  _ __   __ _ _ __(_)_ __   __ _    __ _  __| |(_)_   _ ___| |_ ___  __| | __   _____   _ __ __ ___      __
#  / __/ _ \| '_ ` _ \| '_ \ / _` | '__| | '_ \ / _` |  / _` |/ _` || | | | / __| __/ _ \/ _` | \ \ / / __| | '__/ _` \ \ /\ / /
# | (_| (_) | | | | | | |_) | (_| | |  | | | | | (_| | | (_| | (_| || | |_| \__ \ ||  __/ (_| |  \ V /\__ \ | | | (_| |\ V  V / 
#  \___\___/|_| |_| |_| .__/ \__,_|_|  |_|_| |_|\__, |  \__,_|\__,_|/ |\__,_|___/\__\___|\__,_|   \_/ |___/ |_|  \__,_| \_/\_/  
#                     |_|                       |___/             |__/                                                          
#                                                                                                                      



# LOG-LOG!

# lm for each week
#   week = unique(res$epiweek)[[1]]
r2 = lapply(unique(res$epiweek), function(week) {
  
  tmp = res[res$epiweek == week, ]
  tmp = tmp[!duplicated(tmp$area) , ]
  # fit <- lm(log(true_negative_corrected) ~ log(ne), data = tmp)
  
  tmp$epiweek = paste('Week', tmp$epiweek)
  
  
  
  names(tmp)[names(tmp)=="true_negative"] <- "Raw"
  names(tmp)[names(tmp)=="true_negative_corrected"] <- "TPR-adjusted"
  
  tmp = reshape2::melt(tmp, id.vars = names(tmp)[!names(tmp) %in% c("Raw" , "TPR-adjusted"  )])
  
  
  
  tmp_unadj = tmp[tmp$variable == "Raw", ]
  tmp_adj = tmp[tmp$variable == "TPR-adjusted", ]
  
  
  # sum(abs(summary(lm(log(value)~log(ne), tmp_unadj))$residuals)^2)
  # sum(abs( summary(lm(log(value)~log(ne), tmp_adj))$residuals)^2)
  
  my.formula_unadj <- log(tmp_unadj$value) ~ log(tmp_unadj$ne)
  my.formula_adj <- log(tmp_adj$value) ~ log(tmp_adj$ne)
  
  # red_ss <- sum((mean(log(tmp_unadj$value)) - log(tmp_unadj$value))^2)
  # red_res <-sum(abs(summary(lm(log(value)~log(ne), tmp_unadj))$residuals)^2)
  # blue_ss <- sum((mean(log(tmp_adj$value)) - log(tmp_adj$value))^2)
  # blue_res <-sum(abs(summary(lm(log(value)~log(ne), tmp_adj))$residuals)^2)
  # 
  # 1-red_res/red_ss
  # 1-blue_res/blue_ss
  
  
  pl <- ggplot(tmp, aes(log(ne), log(value), label = area, col = variable, fill = variable)) +  theme_bw()+ geom_point(size = 3)+
    stat_smooth(method = "lm")+
    # geom_label_repel(alpha = 0.2) + geom_errorbarh(aes(xmax = log(neub), xmin = log(nelb))) + 
    # geom_label_repel(data = tmp[tmp$area == "Kent and Medway",], alpha = 0.7, arrow = NULL, nudge_x = log(0.496),
    #                  segment.color = NA) + 
    geom_errorbarh(aes(xmax = log(neub), xmin = log(nelb))) + 
    facet_wrap(~epiweek) + 
    theme(strip.text.x = element_text(size = 20), legend.position = "") + 
    labs(y = "TPR-adjusted SGTF (log scale)", x = "Effective population size (log scale)")+ 
    theme(axis.text=element_text(size=16),
          axis.title=element_text(size=14,face="bold"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank() )   + 
    ggpmisc::stat_poly_eq(data = tmp_unadj, formula = my.formula_unadj, size = 5,
                          mapping = aes(log(ne), log(value),  col = variable, fill = variable,
                                        label = paste( ..adj.rr.label.., sep = "~~~")),  parse = TRUE)+  
    ggpmisc::stat_poly_eq(data = tmp_adj, label.y = 0.85, size = 5,
                          formula = my.formula_adj, mapping = aes(log(ne), log(value), col = variable, fill = variable,
                                                                  label = paste( ..adj.rr.label.., sep = "~~~")),  parse = TRUE)
  
  
  
  
  # fit2 <- lm(log(true_negative_corrected) ~ log(ne), data = tmp)
  
  
  return(list(pl = pl))#,# r2_raw = summary(fit)$adj.r.squared, 
  #r2_corrected = summary(fit2)$adj.r.squared, 
  # p_raw = summary(fit)$coefficients[,4][2], 
  #p_corrected = summary(fit2)$coefficients[,4][2]  ))
}
)




plist = lapply(r2, '[[', 1)
# library(gridExtra)
# n <- length(plist)
# nCol <- floor(sqrt(n))
# do.call("grid.arrange", c(plist, ncol=nCol))


library(grid)
library(gridExtra)


legend <- cowplot::get_legend(
  # create some space to the left of the legend
  plist[[1]] + theme(legend.box.margin = margin(0, 0, 0, 12), legend.position = "top",
                     legend.text = element_text(size=16),
                     legend.title =  element_text(size=20) ) + labs(col = "STGF", fill = "STGF")
)

pl = cowplot::plot_grid(plotlist = plist, align = "v")

pl = cowplot::plot_grid(legend, pl,  align = "v", ncol = 1, rel_heights =  c(0.1, 1))


y.grob <- textGrob("SGTF case numbers (log scale)", 
                   gp=gpar(fontface="bold", fontsize=20), rot=90)

x.grob <- textGrob("Effective population size (log scale)", 
                   gp=gpar(fontface="bold",  fontsize=20))




grid.arrange(arrangeGrob(pl, left = y.grob, bottom = x.grob))
gc()


# save plot 

pdf(file = 'results/SGTF_by_Ne_SG_raw_v_adj.pdf', width = 12, height = 12)
grid.arrange(arrangeGrob(pl, left = y.grob, bottom = x.grob))
dev.off()


png(file = 'results/SGTF_by_Ne_SG_raw_v_adj.png', width = 1000, height = 1000)
grid.arrange(arrangeGrob(pl, left = y.grob, bottom = x.grob))
dev.off()

