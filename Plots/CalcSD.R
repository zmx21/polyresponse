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
resultPath <- '~/bsu_scratch/Random_Forest/rs3821843_rs7340705_rs113210396_rs312487_rs11719824_rs3774530_rs3821856_sbp/'
node_size <- seq(10000,40000,10000)
thresh <- c('1e-5','2e-5','3e-5','4e-5')
comb <- expand.grid(node_size,thresh)
colnames(comb) <- c('node_size','thresh')
n_cores <- 12

nonperm_sd <- pbmclapply(1:nrow(comb),function(i) CalcSD(paste0(resultPath,'0.75_',comb$node_size[i],'_',
                                                                as.character(comb$thresh[i]),'/')),ignore.interactive = T,mc.cores = n_cores)
saveRDS(list(comb = comb,nonperm_sd= nonperm_sd),file=paste0(resultPath,'nonperm_sd.rds'))
