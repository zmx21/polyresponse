source('~/MRC_BSU_Internship/Recursive_Partitioning/InteractionTree.R')
source('~/MRC_BSU_Internship/Recursive_Partitioning/ExtractSubsample.R')
source('~/MRC_BSU_Internship/Recursive_Partitioning/LoadDosage.R')
library(partykit)
library(pbmcapply)
CompareBeta <- function(testing_set,tree){
  #Assign each sample within the testing set to a leaf node
  leafNodeAssignment <- predict(tree,as.data.frame(testing_set$dosageMatrix))
  uniqueLeafNodes <- unique(leafNodeAssignment)
  #Calculate beta coefficients of the subgroups within the testing set.
  testingMainEffects <- rep(NA,length(uniqueLeafNodes))
  names(testingMainEffects) <- uniqueLeafNodes
  for(i in 1:length(uniqueLeafNodes)){
    currentIndicator <- leafNodeAssignment == uniqueLeafNodes[i]
    currentSubset <- lapply(testing_set,function(x) SubsetData(x,currentIndicator))
    testingMainEffects[i] <- FitMainEffectModel(currentSubset$dosageTarget,currentSubset$phenotypes,currentSubset$covariates)$coeff['treatment']
  }
  trainingMainEffects <- lapply(uniqueLeafNodes,function(x) nodeapply(tree,ids = x,FUN = function(n) n$info))
  trainingMainEffects <- sapply(trainingMainEffects,function(x) as.numeric(unlist(strsplit(x = unlist(x),split = "  "))[1]))
  names(trainingMainEffects) <- uniqueLeafNodes
  rse <- sqrt((testingMainEffects-trainingMainEffects)^2)
  names(rse) <- NULL
  return(list(testingMainEffects=testingMainEffects,trainingMainEffects=trainingMainEffects,rse=rse))
}
CalculateBetaError <- function(resultPath,suffix,n_cores){
  data <- readRDS(paste0(resultPath,'data.rds'))
  training_testing_set <- ExtractSubSample(data,readRDS('~/bsu_scratch/UKB_Data/training_set.rds'),readRDS('~/bsu_scratch/UKB_Data/test_set.rds'))
  testing_set <- training_testing_set$outofbag
  rse <- pbmclapply(1:5000,function(x){
    test_tree <- readRDS(paste0(resultPath,suffix,'tree',x,'.rds'))
    CompareBeta(testing_set,test_tree$bootstrapPartyTree)$rse
  },mc.cores = n_cores,ignore.interactive = T)
  return(do.call(c,rse))
}
RunBetaErr <- function(resultPath,suffix,n_cores){
  for(i in 1:length(suffix)){
    betaErr <- CalculateBetaError(resultPath,suffix[i],n_cores)
    saveRDS(betaErr,paste0(resultPath,suffix[i],'beta_err.rds'))
  }
}

n_cores <- 16
args=(commandArgs(TRUE))
node_size <- as.numeric(args[[1]])
resultPath <- '~/bsu_scratch/Random_Forest/rs1262894_sbp/'
RunPredictionFromRF(resultPath,paste0('0.75_',node_size,'_5e-5','/'),n_cores)
RunPredictionFromRF(resultPath,paste0('0.75_',node_size,'_5e-6','/'),n_cores)

resultPath <- '~/bsu_scratch/Random_Forest/rs3821843_rs7340705_rs113210396_rs312487_rs11719824_rs3774530_rs3821856_sbp/'
RunPredictionFromRF(resultPath,paste0('0.75_',node_size,'_5e-5','/'),n_cores)
RunPredictionFromRF(resultPath,paste0('0.75_',node_size,'_5e-6','/'),n_cores)