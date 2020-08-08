#Plot comparison of training set and testing set treatment effects (root squared error)
CalcTrainTestingDiff <- function(path){
  betaErr <- readRDS(paste0(path,'beta_err.rds'))
  return(do.call(c,betaErr))
}

resultPath <- '~/bsu_scratch/LDL_Project_Data/Random_Forest/rs12916_rs17238484_rs5909_rs2303152_rs10066707_rs2006760_LDLdirect/'
node_size <- seq(10000,40000,10000)
thresh <- c('5e-6','1e-5','3e-5','5e-5')
comb <- expand.grid(node_size,thresh)
colnames(comb) <- c('node_size','thresh')

beta_err <- lapply(1:nrow(comb),function(i) CalcTrainTestingDiff(paste0(resultPath,'0.75_',comb$node_size[i],'_',as.character(comb$thresh[i]),'/')))
beta_err_df <- data.frame(p_thresh=numeric(),node_size=numeric(),beta_err=numeric())
for(i in 1:nrow(comb)){
  curBetaErr <- as.numeric(beta_err[[i]])
  beta_err_df <- rbind(beta_err_df,data.frame(p_thresh = rep(as.numeric(as.character(comb$thresh[i])),length(curBetaErr)),
                                      node_size = rep(as.numeric(comb$node_size[i]),length(curBetaErr)),
                                      beta_err = curBetaErr))
}
library(ggplot2)
p <- ggplot(beta_err_df, aes(y=beta_err,x=factor(node_size),fill=factor(node_size))) + 
  geom_boxplot(outlier.shape = NA,position = "dodge") + ylim(0,1.3) + labs(x='Node Size',y='Prediction Error\n(Root Squared Difference between Training and Testing Treatment Effect)') + facet_grid(~factor(p_thresh)) + guides(fill=FALSE)

p <-  geom_bin2d(bins=10)
