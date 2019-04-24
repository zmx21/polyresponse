#mean treatment effects (across all trees in RF). Results saved and used by other plotting functions. 
library(pbmcapply)
CalcMean <- function(path){
  print(path)
  allTrees <- dir(paste0(path,'prediction_betas/'))
  allTrees <- allTrees[sapply(allTrees,function(x) grepl(pattern = 'tree',x=x))]
  
  sumBetas <- rep(0,length(readRDS(paste0(path,'prediction_betas/',allTrees[1]))$nodeAssignment))
  for(i in 1:length(allTrees)){
    curTree <- readRDS(paste0(path,'prediction_betas/',allTrees[i]))
    sumBetas <- sumBetas + curTree$trainingMainEffect[as.character(curTree$nodeAssignment)]
  }
  meanBetas <- sumBetas/(length(allTrees))
  return(meanBetas)
}
resultPath <- '~/bsu_scratch/LDL_Project_Data/Random_Forest/rs12916_rs17238484_rs5909_rs2303152_rs10066707_rs2006760_LDLdirect/'
node_size <- c(1000,5000,10000,20000,30000,40000)
thresh <- c('5e-6','1e-5','3e-5','5e-5','7e-5','1e-4')
comb <- expand.grid(node_size,thresh)
colnames(comb) <- c('node_size','thresh')

mean_betas <- lapply(1:nrow(comb),function(i) CalcMean(paste0(resultPath,'0.75_',comb$node_size[i],'_',as.character(comb$thresh[i]),'/')))
saveRDS(list(comb = comb,mean_betas= mean_betas),file=paste0(resultPath,'mean_betas_training.rds'))
