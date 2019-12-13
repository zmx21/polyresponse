library(ggpubr)
library(pbmcapply)
#hisogram of treatment effects, with different p-value threshold and minimum node sizes.
# mean_betas <- readRDS('~/bsu_scratch/LDL_Project_Data_Aug2019/Random_Forest_Old/rs12916_rs17238484_rs5909_rs2303152_rs10066707_rs2006760_LDLdirect/mean_betas_training.rds')
# mean_betas <- readRDS('~/bsu_scratch/LDL_Project_Data/Random_Forest/rs12916_rs17238484_rs5909_rs2303152_rs10066707_rs2006760_LDLdirect/mea')

n_cores <- 10
node_size <- c(5000,10000,20000,30000,40000)
thresh <- c('7e-5','5e-5','3e-5') #c('6e-6')
comb <- expand.grid(node_size,thresh)
# node_size <- c(5000,node_size)
# thresh <- c('9e-6','7e-6','5e-6','3e-6','9e-5','7e-5','5e-5','3e-5','1e-5')
# comb <- rbind(comb,expand.grid(node_size,thresh))
colnames(comb) <- c('node_size','thresh')

resultPath <- '~/bsu_scratch/LDL_Project_Data_Aug2019/Random_Forest_Old/rs12916_rs17238484_rs5909_rs2303152_rs10066707_rs2006760_LDLdirect/'

# mean_betas_training <- list()
# mean_betas_training$comb <- comb
# mean_betas_training$mean_betas <- list()
# 
# mean_betas_testing <- list()
# mean_betas_testing$comb <- comb
# mean_betas_testing$mean_betas <- list()
n_trees <- 2000

#for(j in 1:nrow(comb)){
#  print(j)
#  cur_comb_result <- lapply(1:n_trees,function(i) readRDS(paste0('~/bsu_scratch/LDL_Project_Data_Aug2019/Random_Forest_Old/rs12916_rs17238484_rs5909_rs2303152_rs10066707_rs2006760_LDLdirect/0.75_',as.character(comb$node_size[j]),'_',as.character(comb$thresh[j]),'_5e-2/prediction_betas/tree',i,'.rds')))

# betas_testing <- lapply(cur_comb_result,function(x) x$testingMainEffects[match(x$nodeAssignment,as.numeric(names(x$testingMainEffects)))])
#  betas_training <- lapply(cur_comb_result,function(x) x$trainingMainEffect[match(x$nodeAssignment,as.numeric(names(x$trainingMainEffect)))])
  # saveRDS(list(cur_comb_result=cur_comb_result,betas_testing=betas_testing,betas_training=betas_training),file = paste0('~/bsu_scratch/LDL_Project_Data_Aug2019/Random_Forest_Old/rs12916_rs17238484_rs5909_rs2303152_rs10066707_rs2006760_LDLdirect/0.75_',as.character(comb$node_size[j]),'_',as.character(comb$thresh[j]),'_5e-2/prediction_betas/betas.rds'))

#  mean_betas_testing <- Reduce('+',betas_testing) / length(betas_testing)
#  mean_betas_training <- Reduce('+',betas_training) / length(betas_testing)
#  saveRDS(list(mean_betas_testing=mean_betas_testing,mean_betas_training=mean_betas_training),file = paste0('~/bsu_scratch/LDL_Project_Data_Aug2019/Random_Forest_Old/rs12916_rs17238484_rs5909_rs2303152_rs10066707_rs2006760_LDLdirect/0.75_',as.character(comb$node_size[j]),'_',as.character(comb$thresh[j]),'_5e-2/prediction_betas/mean_betas.rds'))

#}

for(j in 1:nrow(comb)){
   print(j)
   i = 1
   cur_betas <- readRDS(paste0('~/bsu_scratch/LDL_Project_Data_Aug2019/Random_Forest_Old/rs12916_rs17238484_rs5909_rs2303152_rs10066707_rs2006760_LDLdirect/0.75_',as.character(comb$node_size[j]),'_',as.character(comb$thresh[j]),'_5e-2/prediction_betas_perm/testing_indiv_main_effects_tree',i,'.rds'))
   for(i in 2:n_trees){
     cur_betas <- Map("+",cur_betas,readRDS(paste0('~/bsu_scratch/LDL_Project_Data_Aug2019/Random_Forest_Old/rs12916_rs17238484_rs5909_rs2303152_rs10066707_rs2006760_LDLdirect/0.75_',as.character(comb$node_size[j]),'_',as.character(comb$thresh[j]),'_5e-2/prediction_betas_perm/testing_indiv_main_effects_tree',i,'.rds')))
   }
   mean_betas <- lapply(cur_betas,function(x)x / n_trees)
   saveRDS(mean_betas,file = paste0('~/bsu_scratch/LDL_Project_Data_Aug2019/Random_Forest_Old/rs12916_rs17238484_rs5909_rs2303152_rs10066707_rs2006760_LDLdirect/0.75_',as.character(comb$node_size[j]),'_',as.character(comb$thresh[j]),'_5e-2/prediction_betas_perm/mean_betas.rds'))
   
}

# beta_df <- data.frame(p_thresh=numeric(),node_size=numeric(),beta=numeric())
# for(i in 1:nrow(mean_betas_testing$comb)){
#   curBeta <- as.numeric(mean_betas_testing$mean_betas[[i]])
#   beta_df <- rbind(beta_df,data.frame(p_thresh = rep(as.numeric(as.character(mean_betas_testing$comb$thresh[i])),length(curBeta)),
#                                   node_size = rep(as.numeric(mean_betas_testing$comb$node_size[i]),length(curBeta)),
#                                   beta = curBeta))
# }
# library(ggplot2)
# library(dplyr)
# beta_df <- beta_df %>% dplyr::filter(node_size > 5000 & p_thresh != 3e-5 & p_thresh != 7e-5)
# beta_df$p_thresh <- paste0('p < ',beta_df$p_thresh)
# beta_df$node_size <- paste0('Node Size > ',beta_df$node_size)
# beta_df$p_thresh <- factor(beta_df$p_thresh,levels = unique(beta_df$p_thresh))
# beta_df$node_size <- factor(beta_df$node_size,levels = unique(beta_df$node_size))
# 
# p = ggplot(beta_df, aes(x=beta)) +
#   geom_histogram(bins=50) + ylim(0,40000) + 
#   facet_grid( node_size ~ p_thresh) + labs(x = 'Individualized Treatment Effect', y = 'Count') + theme_pubr(base_size = 14,base_family = 'Myriad Pro',border = T)
# 
# 
# sd_df <- data.frame(p_thresh=numeric(),node_size=numeric(),sd=numeric())
# for(i in 1:nrow(mean_betas$comb)){
#   curBeta <- as.numeric(mean_betas$mean_betas[[i]])
#   sd_df <- rbind(sd_df,data.frame(p_thresh = as.numeric(as.character(mean_betas$comb$thresh[i])),
#                                       node_size = as.numeric(mean_betas$comb$node_size[i]),
#                                       sd = sd(curBeta)))
# }
# sd_df <- sd_df %>% dplyr::filter(node_size > 5000& p_thresh != 3e-5 & p_thresh != 7e-5)
# sd_df$p_thresh <- paste0('p < ',sd_df$p_thresh)
# sd_df$node_size <- paste0('Node Size > ',sd_df$node_size)
# sd_df$p_thresh <- factor(sd_df$p_thresh,levels = unique(sd_df$p_thresh))
# sd_df$node_size <- factor(sd_df$node_size,levels = unique(sd_df$node_size))
# sd_df$sd <- signif(sd_df$sd,3)
# lab = paste0('sd = ',as.character(sapply(sd_df$sd,function(x) formatC(x,digits = 1,format = 'e'))))
# ann_df <- data.frame(p_thresh=sd_df$p_thresh,node_size = sd_df$node_size)
# p + geom_text(data = ann_df,mapping = aes(x=Inf,y=Inf,label = lab),hjust = 1.2,vjust=1.7,size=3.7,family='Myriad Pro') + geom_hline(color = 'white',yintercept = 0,size=1.2)
