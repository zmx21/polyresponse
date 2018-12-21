####################################################################################
#Apply RF to a list of permuted phenotype vector 
#Input: p_value treshold (of predictors to include), and minimum node size of tree.
#Output: .rds files containing predicted treatment effect of each leaf node
####################################################################################
source('~/MRC_BSU_Internship/Recursive_Partitioning/ExtractSubsample.R')
source('~/MRC_BSU_Internship/Recursive_Partitioning/InteractionTree.R')

library(partykit)
library(parallel)
library(pbmcapply)
#Generate training and testing set prediction, given a phenotype matrix
GenerateTestingPrediction <- function(tree,permPhenotypeVect,testingSetSamples){
  nodeAssignment <- predict(tree,as.data.frame(testingSetSamples$dosageMatrix))
  nodeAssignment <- as.vector(nodeAssignment)
  uniqueLeafNodes <- sort(unique(nodeAssignment),decreasing = F)
  #Calculate treatment effect estimates based on testing set
  testingMainEffects <- rep(NA,length(uniqueLeafNodes))
  names(testingMainEffects) <- uniqueLeafNodes
  for(i in 1:length(uniqueLeafNodes)){
    currentIndicator <- nodeAssignment == uniqueLeafNodes[i]
    testingMainEffects[i] <- FitMainEffectModel(testingSetSamples$dosageTarget[currentIndicator],permPhenotypeVect[currentIndicator],testingSetSamples$covariates[currentIndicator,])$coeff['treatment']
  }
  mainEffectPred <- testingMainEffects[as.character(nodeAssignment)]
  names(mainEffectPred) <- NULL
  return(list(mainEffectPred=mainEffectPred,testingMainEffects=testingMainEffects))
}

#Runs predictions of permuted samples in chunks
PermutedPhenoPrediction <- function(testingSetSamples,randomForestPath,n_cores,tree_chunks,permPhenotypeVect){
  system(paste0('mkdir -p ',randomForestPath,'/prediction_betas_pheno_perm/'))
  
  #Parse tree chunk argument
  tree_chunks <- unlist(strsplit(tree_chunks,':'))
  tree_chunks <- tree_chunks[1]:tree_chunks[2]
  
  #Initialize summation matrix
  sumBeta <- lapply(1:length(permPhenotypeVect),function(x) rep(0,nrow(testingSetSamples$dosageMatrix)))
  #Loop through all trees
  pb <- progressBar()
  for(i in 1:(length(tree_chunks))){
    setTxtProgressBar(pb,i/length(tree_chunks))
    #Load trees
    treeObj <- readRDS(paste0(randomForestPath,'tree',tree_chunks[i],'.rds'))
    #Predict treatment effects across bootstrap samples
    permResults <- list()
    if(n_cores == 20){
       chunks <- splitIndices(length(permPhenotypeVect),25)
    }else{
       chunks <- splitIndices(length(permPhenotypeVect),32)
    }
    for(j in 1:length(chunks)){
      permResults <- c(permResults,mclapply(chunks[[j]],function(k){
        GenerateTestingPrediction(treeObj$bootstrapPartyTree,permPhenotypeVect[[k]],testingSetSamples)
      },mc.cores = n_cores))
    }
    sumBeta <- mapply("+",sumBeta,lapply(permResults,function(x) x$mainEffectPred),SIMPLIFY = F)
    testingMainEffects <- lapply(permResults,function(x) x$testingMainEffects)
    saveRDS(testingMainEffects,file=paste0(randomForestPath,'/prediction_betas_pheno_perm/','testing_main_effects_tree',tree_chunks[i],'.rds'))
    remove(permResults);gc(verbose = F);
  }
  close(pb)
  sumBeta <- do.call(rbind,sumBeta)
  saveRDS(sumBeta,file = paste0(randomForestPath,'/prediction_betas_pheno_perm/','sum_beta_',range(tree_chunks)[1],'_',range(tree_chunks)[2],'.rds'))
}

#MAIN FUNCTION
RunPermutedPhenoPrediction <- function(resultPath,suffix,p_thresh,n_cores,tree_chunks){
  #Load training set samples
  data <- readRDS(paste0(resultPath,'data_p_',p_thresh,'.rds'))
  training_testing_set <- ExtractSubSample(data,readRDS('~/bsu_scratch/UKB_Data/training_set.rds'),readRDS('~/bsu_scratch/UKB_Data/test_set.rds'))
  remove(data);gc(verbose = F);
  testing_set <- training_testing_set$outofbag
  remove(training_testing_set);gc(verbose = F);
  #Load permuted phenotype vectors
  permPhenotypeVect <- readRDS(paste0(resultPath,'perm_phenotype_p_',p_thresh,'.rds'))
  PermutedPhenoPrediction(testing_set,paste0(resultPath,suffix),n_cores,tree_chunks,permPhenotypeVect)
}

args=(commandArgs(TRUE))
node_size <- as.numeric(args[[1]])
thresh <- args[[2]]
n_cores <- as.numeric(args[[3]])
tree_chunks <- args[[4]]
# node_size <- 40000
# thresh <- '4e-5'
# n_cores <- 1
# tree_chunks <- '1:1'
print(c('node_size'=node_size,'thresh'=thresh,'n_cores'=n_cores,'tree_chunks'=tree_chunks))
resultPath <- '~/bsu_scratch/Random_Forest/rs3821843_rs7340705_rs113210396_rs312487_rs11719824_rs3774530_rs3821856_sbp/'
RunPermutedPhenoPrediction(resultPath,paste0('0.75_',node_size,'_',thresh,'/'),as.numeric(thresh),n_cores,tree_chunks)
