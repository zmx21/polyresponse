####################################################################################
#Calculates training set treatment effects, and calcualtes the weighted standard deviation
#Input: path to where random forests and dosage data is stored,the thresh used to construct the RF, and 
#whether the true or permuted testing set should be used. 
#Output: list of numerics (weighted cv)
####################################################################################

source('~/MRC_BSU_Internship_LDL/Recursive_Partitioning/InteractionTree.R')
source('~/MRC_BSU_Internship_LDL/Recursive_Partitioning/ExtractSubsample.R')
source('~/MRC_BSU_Internship_LDL/Recursive_Partitioning/LoadDosage.R')
library(partykit)
library(pbmcapply)
#Compare training set and true testing set treatment effects
CompareBeta <- function(testing_set,tree){
  #Assign each sample within the testing set to a leaf node
  leafNodeAssignment <- predict(tree,as.data.frame(testing_set$dosageMatrix))
  uniqueLeafNodes <- sort(unique(leafNodeAssignment),decreasing = F)
  
  #Calculate beta coefficients of the subgroups within the testing set.
  testingBeta <- rep(NA,length(uniqueLeafNodes))
  names(testingBeta) <- uniqueLeafNodes
  testingSampleSizes <- rep(0,length(uniqueLeafNodes))
  names(testingSampleSizes) <- uniqueLeafNodes
  #Loop through all leaves
  for(i in 1:length(uniqueLeafNodes)){
    currentIndicator <- leafNodeAssignment == uniqueLeafNodes[i]
    testingSampleSizes[i] <- sum(currentIndicator)
    testingBeta[i] <- FitMainEffectModel(testing_set$dosageTarget[currentIndicator],testing_set$phenotypes[currentIndicator],testing_set$covariates[currentIndicator,])$coeff['treatment']
  }
  
  #Get training set main effects
  trainingSetInfo <- lapply(uniqueLeafNodes,function(x) nodeapply(tree,ids = x,FUN = function(n) n$info))
  trainingSampleSizes <- as.numeric(sapply(trainingSetInfo,function(x) as.numeric(unlist(strsplit(x = unlist(x),split = "  "))[2])))
  trainingBeta <- as.numeric(sapply(trainingSetInfo,function(x) as.numeric(unlist(strsplit(x = unlist(x),split = "  "))[1])))
  names(trainingBeta) <- uniqueLeafNodes
  names(trainingSampleSizes) <- uniqueLeafNodes
  
  #Calculate weighted cv of training and testing set
  cvTraining <- sqrt(sum(trainingSampleSizes * (trainingBeta - weighted.mean(trainingBeta,trainingSampleSizes))^2)/sum(trainingSampleSizes)) / weighted.mean(trainingBeta,trainingSampleSizes)
  cvTesting <- sqrt(sum(testingSampleSizes * (testingBeta - weighted.mean(testingBeta,testingSampleSizes))^2)/sum(testingSampleSizes)) / weighted.mean(testingBeta,testingSampleSizes)
  return(list(cvTraining=cvTraining,cvTesting=cvTesting))
}
#Compare training set and PERMUTED testing set treatment effects. 
#Precalculated list of permuted predictions should be provided (saved by PermutedPrediction.R)
CompareBetaPerm <- function(perm_results,tree){
  uniqueLeafNodes <- unique(nodeids(tree,terminal=T))
  uniqueLeafNodes <- sort(uniqueLeafNodes,decreasing = F)
  #Get training set main effects
  trainingInfo <- lapply(uniqueLeafNodes,function(x) nodeapply(tree,ids = x,FUN = function(n) n$info))
  trainingSampleSizes <- as.numeric(sapply(trainingInfo,function(x) as.numeric(unlist(strsplit(x = unlist(x),split = "  "))[2])))
  names(trainingSampleSizes) <- uniqueLeafNodes
  cvPerm <- lapply(perm_results,function(x) sqrt(sum(trainingSampleSizes * (x - weighted.mean(x,trainingSampleSizes))^2)/sum(trainingSampleSizes)) / weighted.mean(x,trainingSampleSizes)) 
  return(cvPerm)
}
#Get treatment effects for all trees within random forest and collapse.
CalculateBetaError <- function(resultPath,suffix,p_thresh,n_cores,perm,chunks){
  chunks <- as.numeric(strsplit(chunks,':')[[1]][1]):as.numeric(strsplit(chunks,':')[[1]][2])
  if(perm == 'None'){
    data <- readRDS(paste0(resultPath,'data_p_',p_thresh,'.rds'))
    training_testing_set <- ExtractSubSample(data,readRDS('~/bsu_scratch/LDL_Project_Data/Genotype_Data/training_set.rds'),
                                             readRDS('~/bsu_scratch/LDL_Project_Data/Genotype_Data/test_set.rds'))
    testing_set <- training_testing_set$outofbag
    cv <- lapply(chunks,function(i){
      test_tree <- readRDS(paste0(resultPath,suffix,'tree',i,'.rds'))
      CompareBeta(testing_set,test_tree$bootstrapPartyTree)
    })
    #},mc.cores = n_cores,ignore.interactive = T)
    return(cv)
  }
  else if(perm=='Pheno'){
    cv <- lapply(chunks,function(i){
      #Precalculated list of permuted predictions should be saved (saved by PermutedPrediction.R)
      perm_result <- readRDS(paste0(resultPath,suffix,'prediction_betas_pheno_perm/testing_main_effects_tree',i,'.rds'))
      test_tree <- readRDS(paste0(resultPath,suffix,'tree',i,'.rds'))$bootstrapPartyTree
      CompareBetaPerm(perm_result,test_tree)
    })
      #},mc.cores = n_cores,ignore.interactive = T)
    return(cv)
  }
}

#MAIN FUNCTION. If non-permuted testing set is to be used, calculation is on the fly. 
#If permuted set to be used, precalculated list of permuted predictions should be saved (saved by PermutedPrediction.R)

RunBetaErr <- function(resultPath,suffix,p_thresh,n_cores,perm,chunks='1:2000'){
  file_suffix <- ''
  if(chunks != '1:2000'){
    file_suffix <- paste0('_',strsplit(chunks,':')[[1]][1],'_',strsplit(chunks,':')[[1]][2])
  }
  if(perm=='None'){
    betaErr <- CalculateBetaError(resultPath,suffix,p_thresh,n_cores,perm,chunks)
    saveRDS(betaErr,paste0(resultPath,suffix,'cv',file_suffix,'.rds'))
  }else if (perm =='Pheno'){
    betaErr <- CalculateBetaError(resultPath,suffix,p_thresh,n_cores,perm,chunks)
    saveRDS(betaErr,paste0(resultPath,suffix,'cv_pheno_perm',file_suffix,'.rds'))
  }
}
# resultPath <- '~/bsu_scratch/LDL_Project_Data/Random_Forest/rs12916_rs17238484_rs5909_rs2303152_rs10066707_rs2006760_LDLdirect/'
# args=(commandArgs(TRUE))
# node_size <- as.numeric(args[[1]])
# thresh <- args[[2]]
# n_cores <- as.numeric(args[[3]])
# tree_chunks <- args[[4]]
# perm_type <- args[[5]]
resultPath <- '~/bsu_scratch/LDL_Project_Data/Random_Forest/rs12916_rs17238484_rs5909_rs2303152_rs10066707_rs2006760_LDLdirect/'
node_size=c(5000,10000,20000,30000,40000)
thresh <- c('5e-6','1e-5','3e-5','5e-5','7e-5','1e-4')
comb <- expand.grid(node_size,thresh)
colnames(comb) <- c('node_size','thresh')
n_cores <- 12
for(i in 1:nrow(comb)){
  print(comb[i,])
  RunBetaErr(resultPath,paste0('0.75_',comb$node_size[i],'_',as.character(comb$thresh[i]),'/'),as.numeric(as.character(comb$thresh[i])),n_cores,'None',"1:500")
  RunBetaErr(resultPath,paste0('0.75_',comb$node_size[i],'_',as.character(comb$thresh[i]),'/'),as.numeric(as.character(comb$thresh[i])),n_cores,'None',"501:1000")
  RunBetaErr(resultPath,paste0('0.75_',comb$node_size[i],'_',as.character(comb$thresh[i]),'/'),as.numeric(as.character(comb$thresh[i])),n_cores,'None',"1001:1500")
  RunBetaErr(resultPath,paste0('0.75_',comb$node_size[i],'_',as.character(comb$thresh[i]),'/'),as.numeric(as.character(comb$thresh[i])),n_cores,'None',"1501:2000")
  RunBetaErr(resultPath,paste0('0.75_',comb$node_size[i],'_',as.character(comb$thresh[i]),'/'),as.numeric(as.character(comb$thresh[i])),n_cores,'Pheno',"1:500")
  RunBetaErr(resultPath,paste0('0.75_',comb$node_size[i],'_',as.character(comb$thresh[i]),'/'),as.numeric(as.character(comb$thresh[i])),n_cores,'Pheno',"501:1000")
  RunBetaErr(resultPath,paste0('0.75_',comb$node_size[i],'_',as.character(comb$thresh[i]),'/'),as.numeric(as.character(comb$thresh[i])),n_cores,'Pheno',"1001:1500")
  RunBetaErr(resultPath,paste0('0.75_',comb$node_size[i],'_',as.character(comb$thresh[i]),'/'),as.numeric(as.character(comb$thresh[i])),n_cores,'Pheno',"1501:2000")
}
