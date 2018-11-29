mean_betas <- readRDS('bsu_scratch/Random_Forest/rs3821843_rs7340705_rs113210396_rs312487_rs11719824_rs3774530_rs3821856_sbp/mean_betas.rds')
beta_df <- data.frame(p_thresh=numeric(),node_size=numeric(),beta=numeric())
for(i in 1:nrow(mean_betas$comb)){
  curBeta <- as.numeric(mean_betas$mean_betas[[i]])
  beta_df <- rbind(beta_df,data.frame(p_thresh = rep(as.numeric(as.character(mean_betas$comb$thresh[i])),length(curBeta)),
                                  node_size = rep(as.numeric(mean_betas$comb$node_size[i]),length(curBeta)),
                                  beta = curBeta))
}
library(ggplot2)
p = ggplot(beta_df, aes(x=beta)) +
  geom_histogram() +
  facet_grid( node_size ~ p_thresh) + labs(x = 'Main Treatment Effect')
