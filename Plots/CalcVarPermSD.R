#standard deviation of training set treatment effects. saved and used by other plotting functions. 
library(pbmcapply)
CalcSD <- function(path,ntrees){
  sumBetas <- rep(0,length(readRDS(paste0(path,'prediction_betas/tree1.rds'))$nodeAssignment))
  for(i in 1:ntrees){
    curTree <- readRDS(paste0(path,'prediction_betas/tree',i,'.rds'))
    sumBetas <- sumBetas + curTree$trainingMainEffect[as.character(curTree$nodeAssignment)]
  }
  meanBetas <- sumBetas/ntrees
  sdBetas <- sd(meanBetas)
  return(sdBetas)
}

resultPath <- '~/bsu_scratch/LDL_Project_Data/Random_Forest/Variable_Perm/rs12916_rs17238484_rs5909_rs2303152_rs10066707_rs2006760_LDLdirect/'
node_size <- c(10000,20000,30000,40000)
thresh <- c('5e-6')
comb <- expand.grid(node_size,thresh)
colnames(comb) <- c('node_size','thresh')
n_perm <- 50
ntrees <- 2000

result <- lapply(1:nrow(comb),function(i) unlist(pbmclapply(1:n_perm,function(j) CalcSD(paste0(resultPath,'0.75_',comb$node_size[i],'_',as.character(comb$thresh[i]),'/perm',j,'/'),ntrees),mc.cores=16,ignore.interactive=T)))
saveRDS(list(comb=comb,result=result),file='~/bsu_scratch/LDL_Project_Data/Random_Forest/Variable_Perm/var_perm_sd.rds')
