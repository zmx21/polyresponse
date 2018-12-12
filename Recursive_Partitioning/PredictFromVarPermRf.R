####################################################################################
#Apply RF to testing set, and save predicted treatment effects.
#Input: p_value treshold (of predictors to include), and minimum node size of tree.
#Output: .rds files containing predicted treatment effect of each leaf node
####################################################################################

library(pbmcapply)
library(partykit)


source('~/MRC_BSU_Internship/Recursive_Partitioning/ExtractSubsample.R')
source('~/MRC_BSU_Internship/Recursive_Partitioning/InteractionTree.R')
source('~/MRC_BSU_Internship/Recursive_Partitioning/LoadDosage.R')

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
PredictFromVarPermRF <- function(testingSetSamples,randomForestPath,n_cores){
  system(paste0('mkdir -p ',randomForestPath,'/prediction_betas/'))
  #find all individual trees in random forest
  treePaths <- dir(randomForestPath)
  treePaths <- treePaths[sapply(treePaths,function(x) grepl('tree',x))]
  treePaths <- paste0(randomForestPath,treePaths)
  
  #Extract genotype matrix from testing set
  genotypeMatrix <- as.data.frame(testingSetSamples$dosageMatrix)
  #Calculate beta coeff for each sample for each tree
  trash <- pbmclapply(1:length(treePaths),function(k){
  # trash <- lapply(1:length(treePaths),function(k){
    treeObj <- readRDS(treePaths[k])
    tree <- treeObj$bootstrapPartyTree
    saveRDS(GeneratePrediction(tree,genotypeMatrix,testingSetSamples),paste0(randomForestPath,'/prediction_betas/','tree',k,'.rds'))
    system(paste0('rm ',treePaths[k]))
  #})
  },mc.cores = n_cores,ignore.interactive = T)
}

#MAIN FUNCTION
RunPredictionFromVarPermRF <- function(resultPath,suffix,p_thresh,perm_n,n_cores){
  interaction_path <- '~/parsed_interaction/CACNA1D_sbp.txt'
  targetRS <-  'rs3821843,rs7340705,rs113210396,rs312487,rs11719824,rs3774530,rs3821856'
  phenotype <- 'sbp'
  targetRS <- unlist(strsplit(targetRS,split = ','))
  data <- LoadDosage(as.numeric(p_thresh),interaction_path,phenotype,0.3,0.05,0.5,targetRS)
  data$dosageMatrix <- readRDS(paste0(resultPath,'/dosage_matrix/','dosage_matrix_perm_',as.numeric(perm_n),'_p_',as.numeric(p_thresh),'.rds'))
  training_testing_set <- ExtractSubSample(data,readRDS('~/bsu_scratch/UKB_Data/training_set.rds'),readRDS('~/bsu_scratch/UKB_Data/test_set.rds'))
  testing_set <- training_testing_set$outofbag
  invisible(sapply(paste0(resultPath,suffix),function(x) PredictFromVarPermRF(testing_set,x,n_cores)))
}
args=(commandArgs(TRUE))
node_size <- as.numeric(args[[1]])
thresh <- args[[2]]
n_cores <- as.numeric(args[[3]])
perm_n <- as.numeric(args[[4]])

# node_size <- 10000
# thresh <- '1e-5'
# n_cores <- 16
# perm_n <- 1

print(c("nodesize"=node_size,"thresh"=thresh,"n_cores"=n_cores))
resultPath <- '~/bsu_scratch/Random_Forest/Variable_Perm/rs3821843_rs7340705_rs113210396_rs312487_rs11719824_rs3774530_rs3821856_sbp/'

RunPredictionFromVarPermRF(resultPath,paste0('0.75_',node_size,'_',thresh,'/perm',perm_n,'/'),as.numeric(thresh),perm_n,n_cores)