perm_sd <- readRDS(file = 'bsu_scratch/Random_Forest/rs3821843_rs7340705_rs113210396_rs312487_rs11719824_rs3774530_rs3821856_sbp/perm_sd.rds')
nonperm_sd <- readRDS(file = 'bsu_scratch/Random_Forest/rs3821843_rs7340705_rs113210396_rs312487_rs11719824_rs3774530_rs3821856_sbp/nonperm_sd.rds')
nonperm_sd$comb <- nonperm_sd$comb[-6,]
nonperm_sd$nonperm_sd <- nonperm_sd$nonperm_sd[-6]

sd_df <- data.frame(p_thresh=numeric(),node_size=numeric(),sd=numeric(),nonperm_sd=numeric())
for(i in 1:nrow(perm_sd$comb)){
  curSd <- as.numeric(perm_sd$perm_sd[[i]])
  sd_df <- rbind(sd_df,data.frame(p_thresh = rep(as.numeric(as.character(perm_sd$comb$thresh[i])),length(curSd)),
                      node_size = rep(as.numeric(perm_sd$comb$node_size[i]),length(curSd)),
                      sd = curSd,
                      nonperm_sd = rep(nonperm_sd$nonperm_sd[[i]],length(curSd))))
}
library(ggplot2)
p = ggplot(sd_df, aes(x=sd)) +
  geom_histogram(bins=50) +  xlim(0,0.2) + geom_vline(aes(xintercept=nonperm_sd),data = sd_df,colour = 'red') + 
  facet_grid( node_size ~ p_thresh) + labs(x = 'Standard Deviation of Main Treatment Effect')
