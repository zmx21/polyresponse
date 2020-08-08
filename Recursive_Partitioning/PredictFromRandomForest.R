####################################################################################
#Apply RF to testing set, and save predicted treatment effects.
#Input: p_value treshold (of predictors to include), and minimum node size of tree.
#Output: .rds files containing predicted treatment effect of each leaf node
####################################################################################

library(pbmcapply)
library(partykit)


source('../Recursive_Partitioning/ExtractSubsample.R')
source('../Recursive_Partitioning/InteractionTree.R')
#Generate training and testing set prediction, given a genotype matrix
GeneratePrediction <- function(tree,genotypeMatrix,testingSetSamples){
  nodeAssignment <- predict(tree,genotypeMatrix)
  nodeAssignment <- as.vector(nodeAssignment)
  #Calculate treatment effect estimates based on training set
  trainingMainEffect <- unlist(nodeapply(tree,ids = unique(nodeAssignment),FUN = function(n) as.numeric(unlist(strsplit(x = unlist(n$info),split = "  "))[1])))
  
  #Calculate treatment effect estimates based on testing set
  uniqueLeafNodes <- unique(nodeAssignment)
  testingMainEffects <- rep(NA,length(uniqueLeafNodes))
  names(testingMainEffects) <- uniqueLeafNodes
  for(i in 1:length(uniqueLeafNodes)){
    currentIndicator <- nodeAssignment == uniqueLeafNodes[i]
    currentSubset <- lapply(testingSetSamples,function(x) SubsetData(x,currentIndicator))
    testingMainEffects[i] <- FitMainEffectModel(currentSubset$dosageTarget,currentSubset$phenotypes,currentSubset$covariates)$coeff['treatment']
  }
  results <- list(nodeAssignment = nodeAssignment,testingMainEffects = testingMainEffects,trainingMainEffect=trainingMainEffect)
  return(results)
}

#Runs predictions of testing set samples
PredictFromRF <- function(testingSetSamples,randomForestPath,n_cores){
  system(paste0('mkdir -p ',randomForestPath,'/prediction_betas/'))
  
  #Extract genotype matrix from testing set
  genotypeMatrix <- as.data.frame(testingSetSamples$dosageMatrix)
  #Calculate beta coeff for each sample for each tree
  trash <- pbmclapply(1:2000,function(k){
  #trash <- lapply(1:2000,function(k){
    #print(k)
    treeObj <- readRDS(paste0(randomForestPath,'tree',k,'.rds'))
    tree <- treeObj$bootstrapPartyTree
    saveRDS(GeneratePrediction(tree,genotypeMatrix,testingSetSamples),paste0(randomForestPath,'prediction_betas/','tree',k,'.rds'))
  #})
  },mc.cores = n_cores,ignore.interactive = T)
}

#MAIN FUNCTION
RunPredictionFromRF <- function(resultPath,suffix,p_thresh,n_cores){
  data <- readRDS(paste0(resultPath,'data_p_',p_thresh,'_maf_5e-2.rds'))
  print('Loaded Data')
  training_testing_set <- ExtractSubSample(data,readRDS('~/bsu_scratch/LDL_Project_Data_Aug2019/Genotype_Data/training_set.rds'),
                                           readRDS('~/bsu_scratch/LDL_Project_Data_Aug2019/Genotype_Data/test_set.rds'))
  testing_set <- training_testing_set$outofbag
  PredictFromRF(testing_set,paste0(resultPath,suffix),n_cores)
}
node_size <- c(10000,20000,30000,40000)
thresh <- c('7.5e-6','7.25e-6','6.75e-6','6.5e-6','6e-6')
comb <- expand.grid(node_size,thresh)
node_size <- c(5000,node_size)
thresh <- c('9e-6','7e-6','5e-6','3e-6','9e-5','7e-5','5e-5','3e-5','1e-5')
comb <- rbind(comb,expand.grid(node_size,thresh))
colnames(comb) <- c('node_size','thresh')

print(c("nodesize"=node_size,"thresh"=thresh,"n_cores"=n_cores))
resultPath <- '~/bsu_scratch/LDL_Project_Data_Aug2019/Random_Forest_Old/rs12916_rs17238484_rs5909_rs2303152_rs10066707_rs2006760_LDLdirect/'

for(i in 1:nrow(comb)){
  RunPredictionFromRF(resultPath,paste0('0.75_',comb$node_size[i],'_',comb$thresh[i],'_5e-2/'),as.numeric(as.character(comb$thresh[i])),n_cores)
}

# betas <- lapply(1:2000,function(i) readRDS(paste0('~/bsu_scratch/LDL_Project_Data_Aug2019/Random_Forest_Old/rs12916_rs17238484_rs5909_rs2303152_rs10066707_rs2006760_LDLdirect/0.75_20000_3e-6_5e-2/prediction_betas/tree',i,'.rds')))
# 
# testing_main_effects <- lapply(betas,function(x) x$testingMainEffects[as.character(x$nodeAssignment)])
# training_main_effects <- lapply(betas,function(x) x$trainingMainEffect[as.character(x$nodeAssignment)])
