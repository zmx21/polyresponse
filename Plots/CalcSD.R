#standard deviation of training set treatment effects. saved and used by other plotting functions. 
library(pbmcapply)
CalcSD <- function(path){
  allTrees <- dir(paste0(path,'prediction_betas/'))
  allTrees <- allTrees[sapply(allTrees,function(x) grepl(pattern = 'tree',x=x))]
  
  sumBetas <- rep(0,length(readRDS(paste0(path,'prediction_betas/',allTrees[1]))$nodeAssignment))
  for(i in 1:length(allTrees)){
    curTree <- readRDS(paste0(path,'prediction_betas/',allTrees[i]))
    sumBetas <- sumBetas + curTree$trainingMainEffect[as.character(curTree$nodeAssignment)]
  }
  meanBetas <- sumBetas/(length(allTrees))
  sdBetas <- sd(meanBetas)
  return(sdBetas)
}
resultPath <- '~/bsu_scratch/LDL_Project_Data/Random_Forest/rs12916_rs17238484_rs5909_rs2303152_rs10066707_rs2006760_LDLdirect/'
node_size <- seq(10000,40000,10000)
thresh <- c('5e-6','1e-5','3e-5','5e-5')
comb <- expand.grid(node_size,thresh)
colnames(comb) <- c('node_size','thresh')
n_cores <- 12

nonperm_sd <- pbmclapply(1:nrow(comb),function(i) CalcSD(paste0(resultPath,'0.75_',comb$node_size[i],'_',
                                                                as.character(comb$thresh[i]),'/')),ignore.interactive = T,mc.cores = n_cores)
saveRDS(list(comb = comb,nonperm_sd= nonperm_sd),file=paste0(resultPath,'nonperm_sd.rds'))
