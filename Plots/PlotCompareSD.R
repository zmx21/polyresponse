nonperm_sd <- readRDS(file = '~/bsu_scratch/Random_Forest/rs3821843_rs7340705_rs113210396_rs312487_rs11719824_rs3774530_rs3821856_sbp/nonperm_sd.rds')
perm_sd <- list()
perm_sd$comb <- do.call(rbind,lapply(1:12,function(i) readRDS(file = paste0('~/bsu_scratch/Random_Forest/rs3821843_rs7340705_rs113210396_rs312487_rs11719824_rs3774530_rs3821856_sbp/var_perm_sd_',i,'.rds'))$comb))
perm_sd$perm_sd <- lapply(1:12,function(i) unlist(readRDS(file = paste0('~/bsu_scratch/Random_Forest/rs3821843_rs7340705_rs113210396_rs312487_rs11719824_rs3774530_rs3821856_sbp/var_perm_sd_',i,'.rds'))$result))

nonperm_sd$comb <- nonperm_sd$comb[1:12,]
nonperm_sd$nonperm_sd <- nonperm_sd$nonperm_sd[1:12]


sd_df <- data.frame(p_thresh=numeric(),node_size=numeric(),sd=numeric(),nonperm_sd=numeric(),p_val = numeric())
p_vals <- sapply(1:12,function(i) paste0('p=',sum(unlist(perm_sd$perm_sd[i]) > nonperm_sd$nonperm_sd[[i]]) / 100))

for(i in 1:nrow(perm_sd$comb)){
  curSd <- as.numeric(perm_sd$perm_sd[[i]])
  sd_df <- rbind(sd_df,data.frame(p_thresh = rep(as.numeric(as.character(perm_sd$comb$thresh[i])),length(curSd)),
                      node_size = rep(as.numeric(perm_sd$comb$node_size[i]),length(curSd)),
                      sd = curSd,
                      nonperm_sd = rep(nonperm_sd$nonperm_sd[[i]],length(curSd)),
                      p_val = rep(p_vals[i],length(curSd))))
}
library(ggplot2)

p = ggplot(sd_df, aes(x=sd)) +
  geom_histogram(bins=50) +  xlim(0,0.2) + geom_vline(aes(xintercept=nonperm_sd),data = sd_df,colour = 'red') + 
  geom_text(x=0.15,y=30,aes(label=p_val),data = sd_df) + facet_grid( node_size ~ p_thresh) + labs(x = 'Standard Deviation of Main Treatment Effect') 
