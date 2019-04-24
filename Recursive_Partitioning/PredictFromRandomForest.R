####################################################################################
#Apply RF to testing set, and save predicted treatment effects.
#Input: p_value treshold (of predictors to include), and minimum node size of tree.
#Output: .rds files containing predicted treatment effect of each leaf node
####################################################################################

library(pbmcapply)
library(partykit)


source('~/MRC_BSU_Internship_LDL/Recursive_Partitioning/ExtractSubsample.R')
source('~/MRC_BSU_Internship_LDL/Recursive_Partitioning/InteractionTree.R')
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
  data <- readRDS(paste0(resultPath,'data_p_',p_thresh,'.rds'))
  print('Loaded Data')
  training_testing_set <- ExtractSubSample(data,readRDS('~/bsu_scratch/LDL_Project_Data/Genotype_Data/training_set.rds'),
                                           readRDS('~/bsu_scratch/LDL_Project_Data/Genotype_Data/test_set.rds'))
  testing_set <- training_testing_set$outofbag
  PredictFromRF(testing_set,paste0(resultPath,suffix),n_cores)
}
args=(commandArgs(TRUE))
node_size <- as.numeric(args[[1]])
thresh <- args[[2]]
#n_cores <- as.numeric(args[[3]])

# node_size <- 10000
# thresh <- '7e-5'
n_cores <- 1

print(c("nodesize"=node_size,"thresh"=thresh,"n_cores"=n_cores))
resultPath <- '~/bsu_scratch/LDL_Project_Data/Random_Forest/rs12916_rs17238484_rs5909_rs2303152_rs10066707_rs2006760_LDLdirect/'
#resultPath <- '~/bsu_scratch/LDL_Project_Data/Random_Forest/rs12916_rs72633963_rs55727654_rs17648121_rs2303152_rs62366588_rs75240579_rs111353455_LDLdirect/'

RunPredictionFromRF(resultPath,paste0('0.75_',node_size,'_',thresh,'/'),as.numeric(thresh),n_cores)
