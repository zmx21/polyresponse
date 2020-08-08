#standard deviation of permuted treatment effects. saved and used by other plotting functions. 
library(pbmcapply)
CalcSDPerm <- function(path,node_size,thresh){
  print(path)
  chunks <- list(c(1,500),c(501,1000),c(1001,1500),c(1501,2000))
  allChunks <- sapply(chunks,function(x) paste0(path,'prediction_betas_pheno_perm/sum_beta_',as.character(x[1]),'_',as.character(x[2]),'.rds'))
  
  firstBeta <- readRDS(allChunks[1])
  sumBetas <- firstBeta
  if(length(allChunks) > 1){
    for(i in 2:length(allChunks)){
      print(i)
      curBeta <- readRDS(allChunks[i])
      sumBetas <- sumBetas + curBeta
    }
  }
  meanBeta <- sumBetas / max(unlist(chunks))
  sdBeta <- apply(meanBeta,1,sd)
  return(sdBeta)
}
resultPath <- '~/bsu_scratch/LDL_Project_Data/Random_Forest/rs12916_rs17238484_rs5909_rs2303152_rs10066707_rs2006760_LDLdirect/'
node_size <- c(1000,2500,5000,10000,20000,30000,40000)
thresh <- c('5e-6','1e-5','3e-5','5e-5','7e-5','1e-4')
comb <- expand.grid(node_size,thresh)
colnames(comb) <- c('node_size','thresh')

perm_sd <- lapply(1:nrow(comb),function(i) CalcSDPerm(paste0(resultPath,'0.75_',comb$node_size[i],'_',as.character(comb$thresh[i]),'/'),comb$node_size[i],as.character(comb$thresh[i])))
saveRDS(list(comb = comb,perm_sd= perm_sd),file=paste0(resultPath,'perm_sd.rds'))
