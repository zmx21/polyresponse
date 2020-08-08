####################################################################################
#Calculates and plots individualized treatment effects (testing set)
#Input: p_value treshold (of predictors to include), and minimum node size of tree.
####################################################################################
library(ggpubr)
library(pbmcapply)

n_cores <- 10
node_size <- c(10000,20000,30000,40000)
thresh <- c('3e-6','7e-6','3e-5','7e-5')
comb <- expand.grid(node_size,thresh)
colnames(comb) <- c('node_size','thresh')

resultPath <- '~/bsu_scratch/LDL_Project_Data_Aug2019/Random_Forest_Old/rs12916_rs17238484_rs5909_rs2303152_rs10066707_rs2006760_LDLdirect/'

n_trees <- 2000

for(j in 1:nrow(comb)){
 print(j)
 cur_comb_result <- lapply(1:n_trees,function(i) readRDS(paste0(resultPath,'0.75_',as.character(comb$node_size[j]),'_',as.character(comb$thresh[j]),'_5e-2/prediction_betas/tree',i,'.rds')))

 betas_testing <- lapply(cur_comb_result,function(x) x$testingMainEffects[match(x$nodeAssignment,as.numeric(names(x$testingMainEffects)))])
 betas_training <- lapply(cur_comb_result,function(x) x$trainingMainEffect[match(x$nodeAssignment,as.numeric(names(x$trainingMainEffect)))])
 saveRDS(list(cur_comb_result=cur_comb_result,betas_testing=betas_testing,betas_training=betas_training),file = paste0('~/bsu_scratch/LDL_Project_Data_Aug2019/Random_Forest_Old/rs12916_rs17238484_rs5909_rs2303152_rs10066707_rs2006760_LDLdirect/0.75_',as.character(comb$node_size[j]),'_',as.character(comb$thresh[j]),'_5e-2/prediction_betas/betas.rds'))

 mean_betas_testing <- Reduce('+',betas_testing) / length(betas_testing)
 mean_betas_training <- Reduce('+',betas_training) / length(betas_testing)
 saveRDS(list(mean_betas_testing=mean_betas_testing,mean_betas_training=mean_betas_training),file = paste0(resultPath,'0.75_',as.character(comb$node_size[j]),'_',as.character(comb$thresh[j]),'_5e-2/prediction_betas/mean_betas.rds'))

}

beta_df <- data.frame(p_thresh=numeric(),node_size=numeric(),beta=numeric())
for(i in 1:nrow(comb)){
  curBeta <- readRDS(paste0(resultPath,'0.75_',as.character(comb$node_size[i]),'_',as.character(comb$thresh[i]),'_5e-2/prediction_betas/mean_betas.rds'))$mean_betas_testing
  beta_df <- rbind(beta_df,data.frame(p_thresh = rep(as.numeric(as.character(comb$thresh[i])),length(curBeta)),
                                  node_size = rep(as.numeric(comb$node_size[i]),length(curBeta)),
                                  beta = curBeta))
}
library(ggplot2)
library(dplyr)
beta_df$p_thresh <- paste0('p < ',beta_df$p_thresh)
beta_df$node_size <- paste0('Node Size > ',beta_df$node_size)
beta_df$p_thresh <- factor(beta_df$p_thresh,levels = unique(beta_df$p_thresh))
beta_df$node_size <- factor(beta_df$node_size,levels = unique(beta_df$node_size))

p = ggplot(beta_df, aes(x=beta)) +
  geom_histogram(bins=50) + ylim(0,40000) +
  facet_grid( node_size ~ p_thresh) + labs(x = 'Individual Treatment Effect', y = 'Count')