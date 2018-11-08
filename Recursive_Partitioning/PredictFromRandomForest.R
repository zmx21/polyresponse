source('~/MRC_BSU_Internship/Recursive_Partitioning/ExtractSubsample.R')
library(pbmcapply)
library(partykit)
PredictFromRF <- function(testingSetSamples,randomForestPath,n_cores){
  #find all individual trees in random forest
  treePaths <- dir(randomForestPath)
  treePaths <- treePaths[sapply(treePaths,function(x) grepl('tree',x))]
  treePaths <- paste0(randomForestPath,treePaths)
  
  #Extract genotype matrix from testing set
  genotypeMatrix <- as.data.frame(testingSetSamples$dosageMatrix)
  #Calculate beta coeff for each sample for each tree
  betaCoeff <- pbmclapply(treePaths,function(x){
    treeObj <- readRDS(x)
    tree <- treeObj$bootstrapPartyTree
    nodeids <- predict(tree,genotypeMatrix)
    nodeids <- as.vector(nodeids)
    nodeBetas <- unlist(nodeapply(tree,ids = unique(nodeids),FUN = function(n) as.numeric(unlist(strsplit(x = unlist(n$info),split = "  "))[1])))
    sampleBetas <- sapply(nodeids,function(x) nodeBetas[as.character(x)])
    return(as.vector(sampleBetas))
  },mc.cores = n_cores,ignore.interactive = T)
  betaCoeff <- do.call(rbind,betaCoeff)
  colnames(betaCoeff) <- rownames(genotypeMatrix)
  saveRDS(betaCoeff,paste0(randomForestPath,'rf_predictions.rds'))
}
RunPredictionFromRF <- function(resultPath,suffix,n_cores){
  data <- readRDS(paste0(resultPath,'data.rds'))
  training_testing_set <- ExtractSubSample(data,readRDS('~/bsu_scratch/UKB_Data/training_set.rds'),readRDS('~/bsu_scratch/UKB_Data/test_set.rds'))
  testing_set <- training_testing_set$outofbag
  invisible(sapply(paste0(resultPath,suffix),function(x) PredictFromRF(testing_set,x,n_cores)))
}
n_cores <- 16
args=(commandArgs(TRUE))
node_size <- as.numeric(args[[1]])
print(node_size)
resultPath <- '~/bsu_scratch/Random_Forest/rs1262894_sbp/'
RunPredictionFromRF(resultPath,paste0('0.75_',node_size,'_5e-5','/'),n_cores)
RunPredictionFromRF(resultPath,paste0('0.75_',node_size,'_5e-6','/'),n_cores)

resultPath <- '~/bsu_scratch/Random_Forest/rs3821843_rs7340705_rs113210396_rs312487_rs11719824_rs3774530_rs3821856_sbp/'
RunPredictionFromRF(resultPath,paste0('0.75_',node_size,'_5e-5','/'),n_cores)
RunPredictionFromRF(resultPath,paste0('0.75_',node_size,'_5e-6','/'),n_cores)
