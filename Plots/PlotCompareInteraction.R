CompareInteraction <- function(path){
  true_int <- readRDS(paste0(path,'tree_interactions/interactions.rds'))
  perm_int <- lapply(1:500,function(i) readRDS(paste0(path,'interactions_pheno_perm/interactions_tree',i,'.rds')))
  
  
}


resultPath <- '~/bsu_scratch/LDL_Project_Data/Random_Forest/rs12916_rs17238484_rs5909_rs2303152_rs10066707_rs2006760_LDLdirect/'
node_size <- seq(10000,40000,10000)
thresh <- c('5e-6','1e-5','3e-5','5e-5')
comb <- expand.grid(node_size,thresh)
colnames(comb) <- c('node_size','thresh')
int_diff <- lapply(1:nrow(comb),function(i) CompareInteraction(paste0(resultPath,'0.75_',comb$node_size[i],'_',as.character(comb$thresh[i]),'/')))