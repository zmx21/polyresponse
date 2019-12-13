#mean treatment effects (across all trees in RF). Results saved and used by other plotting functions. 
library(pbmcapply)
CalcMean <- function(path){
  print(path)
  allTrees <- dir(paste0(path,'prediction_betas/'))
  allTrees <- allTrees[sapply(allTrees,function(x) grepl(pattern = 'tree',x=x))]
  
  #Each row is a tree, each column a single sample. Stores predicted treatment effect
  betas <- matrix(NaN,nrow = length(allTrees),ncol = length(readRDS(paste0(path,'prediction_betas/',allTrees[1]))$nodeAssignment))
  for(i in 1:length(allTrees)){
    curTree <- readRDS(paste0(path,'prediction_betas/',allTrees[i]))
    betas[i,] <- curTree$trainingMainEffect[as.character(curTree$nodeAssignment)]
  }
  saveRDS(betas,paste0(path,'pred_beta_mat.rds'))
  return(mean_betas = as.vector(colMeans(betas)))
}
resultPath <- '~/bsu_scratch/LDL_Project_Data_Aug2019/Random_Forest_Old/rs12916_rs17238484_rs5909_rs2303152_rs10066707_rs2006760_LDLdirect/'
node_size <- c(10000,20000,30000,40000)
thresh <- c('7.5e-6','7.25e-6','6.75e-6','6.5e-6','6e-6')
comb <- expand.grid(node_size,thresh)
node_size <- c(5000,node_size)
thresh <- c('9e-6','7e-6','5e-6','3e-6','9e-5','7e-5','5e-5','3e-5','1e-5')
comb <- rbind(comb,expand.grid(node_size,thresh))
colnames(comb) <- c('node_size','thresh')

mean_betas <- lapply(1:nrow(comb),function(i) CalcMean(paste0(resultPath,'0.75_',comb$node_size[i],'_',as.character(comb$thresh[i]),'_5e-2/')))
saveRDS(list(comb = comb,mean_betas= mean_betas),file=paste0(resultPath,'mean_betas_training.rds'))

