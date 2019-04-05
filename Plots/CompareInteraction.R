CompareInteractions <- function(targetRS,suffix,n_perm){
  print(suffix)
  #Get sum of interactions from true tree
  trueTreeInteractions <- readRDS(paste0('~/bsu_scratch/LDL_Project_Data/Random_Forest/',targetRS,suffix,'tree_interactions/interactions.rds'))
  trueMeanSumOfInteractions <- mean(sapply(trueTreeInteractions,function(x) max(x)))
  
  #Get sum of interactions from tree with random predictors
  permTreeInteractions <- lapply(1:n_perm,function(i) readRDS(paste0('~/bsu_scratch/LDL_Project_Data/Random_Forest/Variable_Perm_Interaction/',
                                                             targetRS,suffix,'perm',i,'/interactions.rds')))
  permMeanSumOfInteractions <- sapply(permTreeInteractions,function(x) mean(sapply(x,function(y) max(y))))
  return(list(true=trueMeanSumOfInteractions,perm=permMeanSumOfInteractions))
}


targetRS <- 'rs12916_rs17238484_rs5909_rs2303152_rs10066707_rs2006760_LDLdirect/'
node_size <- seq(10000,40000,10000)
#thresh <- c('5e-6','1e-5','3e-5','5e-5')
thresh <- c('5e-6','1e-5','3e-5')
comb <- expand.grid(node_size,thresh)
colnames(comb) <- c('node_size','thresh')
n_perm <- 50

res <- lapply(1:nrow(comb),function(i) CompareInteractions(targetRS,paste0('0.75_',comb$node_size[i],'_',as.character(comb$thresh[i]),'/'),n_perm))
saveRDS(res,file='~/bsu_scratch/LDL_Project_Data/Random_Forest/Variable_Perm_Interaction/interaction_max_comparison.rds')