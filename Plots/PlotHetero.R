#hisogram of treatment effects, with different p-value threshold and minimum node sizes.
mean_betas <- readRDS('~/bsu_scratch/LDL_Project_Data/Random_Forest/rs12916_rs17238484_rs5909_rs2303152_rs10066707_rs2006760_LDLdirect/mean_betas.rds')
beta_df <- data.frame(p_thresh=numeric(),node_size=numeric(),beta=numeric())
for(i in 1:nrow(mean_betas$comb)){
  curBeta <- as.numeric(mean_betas$mean_betas[[i]])
  beta_df <- rbind(beta_df,data.frame(p_thresh = rep(as.numeric(as.character(mean_betas$comb$thresh[i])),length(curBeta)),
                                  node_size = rep(as.numeric(mean_betas$comb$node_size[i]),length(curBeta)),
                                  beta = curBeta))
}
library(ggplot2)
beta_df$p_thresh <- paste0('p < ',beta_df$p_thresh)
beta_df$node_size <- paste0('Node Size > ',beta_df$node_size)
p = ggplot(beta_df, aes(x=beta)) +
  geom_histogram() +
  facet_grid( node_size ~ p_thresh) + labs(x = 'Individualized Treatment Effect', y = 'Count') + theme(text = element_text(size=14))
