#standard deviation of training set treatment effects. saved and used by other plotting functions. 
library(pbmcapply)
CalcSD <- function(path){
  allTrees <- dir(paste0(path,'prediction_betas/'))
  allTrees <- allTrees[sapply(allTrees,function(x) grepl(pattern = 'tree',x=x))]
  
  sumBetas <- rep(0,length(readRDS(paste0(path,'prediction_betas/',allTrees[1]))$nodeAssignment))
  for(i in 1:length(allTrees)){
    print(i)
    curTree <- readRDS(paste0(path,'prediction_betas/',allTrees[i]))
    sumBetas <- sumBetas + curTree$trainingMainEffect[as.character(curTree$nodeAssignment)]
  }
  meanBetas <- sumBetas/(length(allTrees))
  sdBetas <- sd(meanBetas)
  return(sdBetas)
}

resultPath <- '~/bsu_scratch/Random_Forest/Variable_Perm/rs3821843_rs7340705_rs113210396_rs312487_rs11719824_rs3774530_rs3821856_sbp/'
node_size <- seq(10000,40000,10000)
thresh <- c('1e-5','2e-5','3e-5','4e-5')
comb <- expand.grid(node_size,thresh)
colnames(comb) <- c('node_size','thresh')

comb <- comb[1:2,]
n_perm <- 200


result <- lapply(1:nrow(comb),function(i) unlist(pbmclapply(1:n_perm,function(j) CalcSD(paste0(resultPath,'0.75_',comb$node_size[i],'_',as.character(comb$thresh[i]),'/perm',j,'/')),mc.cores=32,ignore.interactive=T)))

saveRDS(list(comb=comb,result=result),file='~/bsu_scratch/Random_Forest/var_perm_sd.rds')
