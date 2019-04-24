#hisogram of treatment effects, with different p-value threshold and minimum node sizes.
#mean_betas <- readRDS('~/bsu_scratch/LDL_Project_Data/Random_Forest/rs12916_rs17238484_rs5909_rs2303152_rs10066707_rs2006760_LDLdirect/mean_betas.rds')
mean_betas <- readRDS('/media/sf_D_DRIVE/bsu_scratch/Random_Forest/rs12916_rs17238484_rs5909_rs2303152_rs10066707_rs2006760_LDLdirect/mean_betas_testing.rds')

beta_df <- data.frame(p_thresh=numeric(),node_size=numeric(),beta=numeric())
for(i in 1:nrow(mean_betas$comb)){
  curBeta <- as.numeric(mean_betas$mean_betas[[i]])
  beta_df <- rbind(beta_df,data.frame(p_thresh = rep(as.numeric(as.character(mean_betas$comb$thresh[i])),length(curBeta)),
                                  node_size = rep(as.numeric(mean_betas$comb$node_size[i]),length(curBeta)),
                                  beta = curBeta))
}
library(ggplot2)
library(dplyr)
beta_df <- beta_df %>% dplyr::filter(node_size > 5000 & p_thresh != 3e-5 & p_thresh != 7e-5)
beta_df$p_thresh <- paste0('p < ',beta_df$p_thresh)
beta_df$node_size <- paste0('Node Size > ',beta_df$node_size)
beta_df$p_thresh <- factor(beta_df$p_thresh,levels = unique(beta_df$p_thresh))
beta_df$node_size <- factor(beta_df$node_size,levels = unique(beta_df$node_size))

p = ggplot(beta_df, aes(x=beta)) +
  geom_histogram(bins=50) + ylim(0,40000) + 
  facet_grid( node_size ~ p_thresh) + labs(x = 'Individualized Treatment Effect', y = 'Count') + 
  theme(text = element_text(size=14,family = 'Myriad Pro'))


sd_df <- data.frame(p_thresh=numeric(),node_size=numeric(),sd=numeric())
for(i in 1:nrow(mean_betas$comb)){
  curBeta <- as.numeric(mean_betas$mean_betas[[i]])
  sd_df <- rbind(sd_df,data.frame(p_thresh = as.numeric(as.character(mean_betas$comb$thresh[i])),
                                      node_size = as.numeric(mean_betas$comb$node_size[i]),
                                      sd = sd(curBeta)))
}
sd_df <- sd_df %>% dplyr::filter(node_size > 5000& p_thresh != 3e-5 & p_thresh != 7e-5)
sd_df$p_thresh <- paste0('p < ',sd_df$p_thresh)
sd_df$node_size <- paste0('Node Size > ',sd_df$node_size)
sd_df$p_thresh <- factor(sd_df$p_thresh,levels = unique(sd_df$p_thresh))
sd_df$node_size <- factor(sd_df$node_size,levels = unique(sd_df$node_size))
sd_df$sd <- signif(sd_df$sd,3)
lab = paste0('sd = ',as.character(sapply(sd_df$sd,function(x) formatC(x,digits = 1,format = 'e'))))
ann_df <- data.frame(p_thresh=sd_df$p_thresh,node_size = sd_df$node_size)
p + geom_text(data = ann_df,mapping = aes(x=Inf,y=Inf,label = lab),hjust = 1,vjust=1.7,size=3.5,family='Myriad Pro') + geom_hline(color = 'white',yintercept = 0,size=1.2)
