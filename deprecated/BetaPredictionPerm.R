library(pbmcapply)
CalcSD <- function(path){
  chunks <- list(c(1,500),c(501,1000),c(1001,1500),c(1501,2000),c(2001,2500),c(2501,3000),c(3001,4500),c(4501,5000))
  allChunks <- sapply(chunks,function(x) paste0(path,'prediction_betas_perm/sum_beta_',as.character(x[1]),'_',as.character(x[2]),'.rds'))
  
  firstBeta <- readRDS(allChunks[1])
  sumBetas <- firstBeta
  if(length(allChunks) > 1){
    for(i in 2:length(allChunks)){
      curBeta <- readRDS(allChunks[i])
      sumBetas <- sumBetas + curBeta
    }
  }
  meanBeta <- sumBetas / max(unlist(chunks))
  sdBeta <- apply(meanBeta,1,sd)
  return(sdBeta)
}
resultPath <- '~/bsu_scratch/Random_Forest/rs3821843_rs7340705_rs113210396_rs312487_rs11719824_rs3774530_rs3821856_sbp/'
node_size <- seq(10000,40000,10000)
thresh <- c('1e-5','2e-5','3e-5','4e-5')
comb <- expand.grid(node_size,thresh)
colnames(comb) <- c('node_size','thresh')

CalcSD(paste0(resultPath,'0.75_',comb$node_size[1],'_',comb$thresh[1],'/'))
# pbmclapply(1:nrow(comb),function(i) CalcSD(paste0(resultPath,'0.75_',comb$node_size[i],'_',comb$thresh[i],'/')))



