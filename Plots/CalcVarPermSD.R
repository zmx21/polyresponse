#standard deviation of training set treatment effects. saved and used by other plotting functions. 
library(pbmcapply)
CalcSD <- function(path){
  sumBetas <- rep(0,length(readRDS(paste0(path,'prediction_betas/tree1.rds'))$nodeAssignment))
  for(i in 1:5000){
    curTree <- readRDS(paste0(path,'prediction_betas/tree',i,'.rds'))
    sumBetas <- sumBetas + curTree$trainingMainEffect[as.character(curTree$nodeAssignment)]
  }
  meanBetas <- sumBetas/5000
  sdBetas <- sd(meanBetas)
  return(sdBetas)
}
args <- (commandArgs(TRUE))
index <- args[[1]] 
resultPath <- '~/bsu_scratch/Random_Forest/Variable_Perm/rs3821843_rs7340705_rs113210396_rs312487_rs11719824_rs3774530_rs3821856_sbp/'
node_size <- c(10000,20000,30000,40000)
thresh <- c('1e-5','2e-5','3e-5')
#node_size = 20000
#thresh = '1e-5'
comb <- expand.grid(node_size,thresh)
colnames(comb) <- c('node_size','thresh')
n_perm <- 100
comb <- comb[as.numeric(index),]

result <- lapply(1:nrow(comb),function(i) unlist(pbmclapply(1:n_perm,function(j) CalcSD(paste0(resultPath,'0.75_',comb$node_size[i],'_',as.character(comb$thresh[i]),'/perm',j,'/')),mc.cores=26,ignore.interactive=T)))
saveRDS(list(comb=comb,result=result),file=paste0('~/bsu_scratch/Random_Forest/var_perm_sd_',index,'.rds'))
