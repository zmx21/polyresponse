nonperm_sd <- readRDS(file = '~/bsu_scratch/LDL_Project_Data/Random_Forest/rs12916_rs17238484_rs5909_rs2303152_rs10066707_rs2006760_LDLdirect/nonperm_sd.rds')
perm_sd <- readRDS(file='~/bsu_scratch/LDL_Project_Data/Random_Forest/rs12916_rs17238484_rs5909_rs2303152_rs10066707_rs2006760_LDLdirect/perm_sd.rds')

sd_df <- data.frame(p_thresh=numeric(),node_size=numeric(),sd=numeric(),nonperm_sd=numeric(),p_val = numeric())

p_vals <- sapply(1:length(perm_sd$perm_sd),function(i) paste0('p=',sum(abs(unlist(perm_sd$perm_sd[i])) > abs(nonperm_sd$nonperm_sd[[i]])) / 1000))

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
