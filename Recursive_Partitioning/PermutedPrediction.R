source('~/MRC_BSU_Internship/Recursive_Partitioning/ExtractSubsample.R')
source('~/MRC_BSU_Internship/Recursive_Partitioning/InteractionTree.R')

library(partykit)
library(parallel)
#Generate training and testing set prediction, given a genotype matrix
GenerateTestingPrediction <- function(tree,genotypeMatrix,testingSetSamples){
  nodeAssignment <- predict(tree,genotypeMatrix)
  nodeAssignment <- as.vector(nodeAssignment)
  uniqueLeafNodes <- sort(unique(nodeAssignment),decreasing = F)
  #Calculate treatment effect estimates based on testing set
  testingMainEffects <- rep(NA,length(uniqueLeafNodes))
  names(testingMainEffects) <- uniqueLeafNodes
  for(i in 1:length(uniqueLeafNodes)){
    currentIndicator <- nodeAssignment == uniqueLeafNodes[i]
    testingMainEffects[i] <- FitMainEffectModel(testingSetSamples$dosageTarget[currentIndicator],testingSetSamples$phenotypes[currentIndicator],testingSetSamples$covariates[currentIndicator,])$coeff['treatment']
  }
  results <- list(nodeAssignment = nodeAssignment,testingMainEffects = testingMainEffects)
  return(results)
}


PermutedPrediction <- function(testingSetSamples,randomForestPath,n_cores,tree_chunks,permGenotypeMatrix){
  system(paste0('mkdir -p ',randomForestPath,'/prediction_betas_perm/'))
  
  #find all individual trees in random forest
  treePaths <- dir(randomForestPath)
  treePaths <- treePaths[sapply(treePaths,function(x) grepl('tree',x))]
  treePaths <- paste0(randomForestPath,treePaths)
  
  #For each tree, predict using the permuted matrix 
  tree_chunks <- unlist(strsplit(tree_chunks,':'))
  tree_chunks <- tree_chunks[1]:tree_chunks[2]
  tree_chunks <- split(tree_chunks,seq(length(tree_chunks)-1)%/%n_cores)
  
  #Seperate into chunks of number of cores (reduce memory usage)
  for(i in 1:(length(tree_chunks))){
    currentTreeChunk <- tree_chunks[[i]]
    print(currentTreeChunk)
    
    trash <- mclapply(currentTreeChunk,function(k){
      # k <- 1
      treeObj <- readRDS(treePaths[k])
      tree <- treeObj$bootstrapPartyTree
      permResults <- lapply(permGenotypeMatrix,function(x) GenerateTestingPrediction(tree,x,testingSetSamples))
      rm(treeObj);gc(verbose = F);
      curTreeResults <- lapply(c('nodeAssignment','testingMainEffects'),function(i) do.call(rbind,lapply(permResults,function(x) x[[i]])))
      rm(permResults);gc(verbose = F);
      names(curTreeResults) <- c('nodeAssignment','testingMainEffects')
      saveRDS(curTreeResults,paste0(randomForestPath,'/prediction_betas_perm/','tree',k,'.rds'))
      rm(curTreeResults); gc(verbose = F);
      return(NA)
    },mc.preschedule = F,mc.cores = n_cores)
    remove(trash);gc(verbose = T);
  }
}
RunPermutedPrediction <- function(resultPath,suffix,p_thresh,n_cores,tree_chunks){
  #Load training set samples
  data <- readRDS(paste0(resultPath,'data_p_',p_thresh,'.rds'))
  training_testing_set <- ExtractSubSample(data,readRDS('~/bsu_scratch/UKB_Data/training_set.rds'),readRDS('~/bsu_scratch/UKB_Data/test_set.rds'))
  remove(data);gc(verbose = F);
  testing_set <- training_testing_set$outofbag
  remove(training_testing_set);gc(verbose = F);
  #Load permuted genotype matrix
  permGenotypeMatrix <- readRDS(paste0(resultPath,'perm_genotype_p_',p_thresh,'.rds'))
  PermutedPrediction(testing_set,paste0(resultPath,suffix),n_cores,tree_chunks,permGenotypeMatrix)
}

args=(commandArgs(TRUE))
node_size <- as.numeric(args[[1]])
thresh <- args[[2]]
n_cores <- as.numeric(args[[3]])
tree_chunks <- args[[4]]
# node_size <- 10000
# thresh <- '3e-5'
# n_cores <- 16
# tree_chunks <- '1:33'
print(c('node_size'=node_size,'thresh'=thresh,'n_cores'=n_cores,'tree_chunks'=tree_chunks))
resultPath <- '~/bsu_scratch/Random_Forest/rs3821843_rs7340705_rs113210396_rs312487_rs11719824_rs3774530_rs3821856_sbp/'
RunPermutedPrediction(resultPath,paste0('0.75_',node_size,'_',thresh,'/'),as.numeric(thresh),n_cores,tree_chunks)
