SaveTree = function(cur_dir,index){
  cur_obj = readRDS(paste0(cur_dir,'tree',index,'.rds'))
  cur_tree = cur_obj$bootstrapPartyTree
  
  #Set empty dataframe (to remove UKB data)
  cur_tree$data = cur_tree$data[0,]
  cur_tree$fitted = NULL
  
  saveRDS(cur_tree,paste0(cur_dir,'pub_trees/','tree',index,'.rds'))
}

resultPath <- '~/bsu_scratch/LDL_Project_Data_Aug2019/Random_Forest_Old/rs12916_rs17238484_rs5909_rs2303152_rs10066707_rs2006760_LDLdirect/'
node_size <- c(10000,20000,30000,40000)
thresh <- c('7.5e-6','7.25e-6','6.75e-6','6.5e-6','6e-6')
comb <- expand.grid(node_size,thresh)
node_size <- c(5000,node_size)
thresh <- c('9e-6','7e-6','5e-6','3e-6','9e-5','7e-5','5e-5','3e-5','1e-5')
comb <- rbind(comb,expand.grid(node_size,thresh))
colnames(comb) <- c('node_size','thresh')

for(i in 1:nrow(comb)){
  print(i)
  cur_dir = paste0(resultPath,'0.75_',comb$node_size[i],'_',as.character(comb$thresh[i]),'_5e-2/')
  system(paste0('mkdir -p ',cur_dir,'pub_trees'))
  lapply(1:2000,function(i) SaveTree(cur_dir,i))
}
