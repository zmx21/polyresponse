####################################################################################
#Apply RF to a list of permuted genotype matrices (should be saved already by GeneratePermutedGenoMatrix.R)
#Input: p_value treshold (of predictors to include), and minimum node size of tree.
#Output: .rds files containing predicted treatment effect of each leaf node
####################################################################################
source('~/MRC_BSU_Internship/Recursive_Partitioning/ExtractSubsample.R')
source('~/MRC_BSU_Internship/Recursive_Partitioning/InteractionTree.R')

library(partykit)
library(parallel)
library(pbmcapply)
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
  mainEffectPred <- testingMainEffects[as.character(nodeAssignment)]
  names(mainEffectPred) <- NULL
  return(list(mainEffectPred=mainEffectPred,testingMainEffects=testingMainEffects))
}

#Runs predictions of permuted samples in chunks
PermutedPrediction <- function(testingSetSamples,randomForestPath,n_cores,tree_chunks,permGenotypeMatrix){
  system(paste0('mkdir -p ',randomForestPath,'/prediction_betas_perm/'))
  
  #find all individual trees in random forest
  treePaths <- dir(randomForestPath)
  treePaths <- treePaths[sapply(treePaths,function(x) grepl('tree',x))]
  treePaths <- paste0(randomForestPath,treePaths)
  
  #For each tree, predict using the permuted matrix 
  tree_chunks <- unlist(strsplit(tree_chunks,':'))
  tree_chunks <- tree_chunks[1]:tree_chunks[2]

  #Initialize summation matrix
  sumBeta <- lapply(1:length(permGenotypeMatrix),function(x) rep(0,nrow(testingSetSamples$dosageMatrix)))
  #Loop through all trees
  pb <- progressBar()
  for(i in 1:(length(tree_chunks))){
      setTxtProgressBar(pb,i/length(tree_chunks))
      #Load trees
      treeObj <- readRDS(treePaths[tree_chunks[i]])
      #Predict treatment effects across bootstrap samples
      permResults <- list()
      if(n_cores==16){
        chunks <- splitIndices(length(permGenotypeMatrix),10)
      }else{
        chunks <- splitIndices(length(permGenotypeMatrix),50)
      }
      for(j in 1:length(chunks)){
        permResults <- c(permResults,mclapply(chunks[[j]],function(k){
          GenerateTestingPrediction(treeObj$bootstrapPartyTree,permGenotypeMatrix[[k]],testingSetSamples)
        },mc.cores =  n_cores))
      }
      sumBeta <- mapply("+",sumBeta,lapply(permResults,function(x) x$mainEffectPred),SIMPLIFY = F)
      testingMainEffects <- lapply(permResults,function(x) x$testingMainEffects)
      saveRDS(testingMainEffects,file=paste0(randomForestPath,'/prediction_betas_perm/','testing_main_effects_tree',tree_chunks[i],'.rds'))
      remove(permResults);gc(verbose = F);
  }
  close(pb)
  sumBeta <- do.call(rbind,sumBeta)
  saveRDS(sumBeta,file = paste0(randomForestPath,'/prediction_betas_perm/','sum_beta_',range(tree_chunks)[1],'_',range(tree_chunks)[2],'.rds'))
}

#MAIN FUNCTION
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
# node_size <- 40000
# thresh <- '1e-5'
# n_cores <- 16
# tree_chunks <- '2:2'
print(c('node_size'=node_size,'thresh'=thresh,'n_cores'=n_cores,'tree_chunks'=tree_chunks))
resultPath <- '~/bsu_scratch/Random_Forest/rs3821843_rs7340705_rs113210396_rs312487_rs11719824_rs3774530_rs3821856_sbp/'
RunPermutedPrediction(resultPath,paste0('0.75_',node_size,'_',thresh,'/'),as.numeric(thresh),n_cores,tree_chunks)
