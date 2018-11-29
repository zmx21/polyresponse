####################################################################################
#Calculates training set treatment effects, and calcualtes the Root Squared error (against training set)
#Input: path to where random forests and dosage data is stored,the thresh used to construct the RF, and 
#whether the true or permuted testing set should be used. 
#Output: list of numerics (root squared error)
####################################################################################

source('~/MRC_BSU_Internship/Recursive_Partitioning/InteractionTree.R')
source('~/MRC_BSU_Internship/Recursive_Partitioning/ExtractSubsample.R')
source('~/MRC_BSU_Internship/Recursive_Partitioning/LoadDosage.R')
library(partykit)
library(pbmcapply)
#Compare training set and true testing set treatment effects
CompareBeta <- function(testing_set,tree){
  #Assign each sample within the testing set to a leaf node
  leafNodeAssignment <- predict(tree,as.data.frame(testing_set$dosageMatrix))
  uniqueLeafNodes <- sort(unique(leafNodeAssignment),decreasing = F)
  
  #Calculate beta coefficients of the subgroups within the testing set.
  testingMainEffects <- rep(NA,length(uniqueLeafNodes))
  names(testingMainEffects) <- uniqueLeafNodes
  #Loop through all leaves
  for(i in 1:length(uniqueLeafNodes)){
    currentIndicator <- leafNodeAssignment == uniqueLeafNodes[i]
    testingMainEffects[i] <- FitMainEffectModel(testing_set$dosageTarget[currentIndicator],testing_set$phenotypes[currentIndicator],testing_set$covariates[currentIndicator,])$coeff['treatment']
  }
  
  #Get training set main effects
  trainingMainEffects <- lapply(uniqueLeafNodes,function(x) nodeapply(tree,ids = x,FUN = function(n) n$info))
  trainingMainEffects <- sapply(trainingMainEffects,function(x) as.numeric(unlist(strsplit(x = unlist(x),split = "  "))[1]))
  names(trainingMainEffects) <- uniqueLeafNodes
  #Calculate RSE
  rse <- sqrt((testingMainEffects-trainingMainEffects)^2)
  names(rse) <- uniqueLeafNodes
  return(list(rse=rse,trainingMainEffects=trainingMainEffects,testingMainEffects=testingMainEffects))
}
#Compare training set and PERMUTED testing set treatment effects. 
#Precalculated list of permuted predictions should be provided (saved by PermutedPrediction.R)
CompareBetaPerm <- function(perm_results,tree){
  uniqueLeafNodes <- unique(nodeids(tree,terminal=T))
  uniqueLeafNodes <- sort(uniqueLeafNodes,decreasing = F)
  #Get training set main effects
  trainingMainEffects <- lapply(uniqueLeafNodes,function(x) nodeapply(tree,ids = x,FUN = function(n) n$info))
  trainingMainEffects <- sapply(trainingMainEffects,function(x) as.numeric(unlist(strsplit(x = unlist(x),split = "  "))[1]))
  names(trainingMainEffects) <- uniqueLeafNodes
  rse <- lapply(perm_results,function(x) sqrt((x-trainingMainEffects)^2))
  return(rse)
}
#Get treatment effects for all trees within random forest and collapse.
CalculateBetaError <- function(resultPath,suffix,p_thresh,n_cores,perm){
  if(perm){
    rse <- pbmclapply(1:5000,function(i){
      #Perm ordering is different
      #find all individual trees in random forest
      treePaths <- dir(paste0(resultPath,suffix))
      treePaths <- treePaths[sapply(treePaths,function(x) grepl('tree',x))]
      
      #Precalculated list of permuted predictions should be saved (saved by PermutedPrediction.R)
      perm_result <- readRDS(paste0(resultPath,suffix,'prediction_betas_perm/testing_main_effects_tree',which(treePaths==paste0('tree',i,'.rds')),'.rds'))
      test_tree <- readRDS(paste0(resultPath,suffix,'tree',i,'.rds'))$bootstrapPartyTree
      CompareBetaPerm(perm_result,test_tree)
    },mc.cores = n_cores,ignore.interactive = T)
    return(rse)
  }else{
    data <- readRDS(paste0(resultPath,'data_p_',p_thresh,'.rds'))
    training_testing_set <- ExtractSubSample(data,readRDS('~/bsu_scratch/UKB_Data/training_set.rds'),readRDS('~/bsu_scratch/UKB_Data/test_set.rds'))
    testing_set <- training_testing_set$outofbag
    rse <- pbmclapply(1:5000,function(i){
      test_tree <- readRDS(paste0(resultPath,suffix,'tree',i,'.rds'))
      CompareBeta(testing_set,test_tree$bootstrapPartyTree)$rse
    },mc.cores = n_cores,ignore.interactive = T)
    return(rse)
  }
}

#MAIN FUNCTION. If non-permuted testing set is to be used, calculation is on the fly. 
#If permuted set to be used, precalculated list of permuted predictions should be saved (saved by PermutedPrediction.R)

RunBetaErr <- function(resultPath,suffix,p_thresh,n_cores,perm){
  if(!perm){
    betaErr <- CalculateBetaError(resultPath,suffix,p_thresh,n_cores,perm)
    saveRDS(betaErr,paste0(resultPath,suffix,'beta_err.rds'))
  }else{
    betaErr <- CalculateBetaError(resultPath,suffix,p_thresh,n_cores,perm)
    saveRDS(betaErr,paste0(resultPath,suffix,'beta_err_perm.rds'))
  }
}
resultPath <- '~/bsu_scratch/Random_Forest/rs3821843_rs7340705_rs113210396_rs312487_rs11719824_rs3774530_rs3821856_sbp/'
node_size <- seq(10000,40000,10000)
thresh <- c('1e-5','2e-5','3e-5','4e-5')
comb <- expand.grid(node_size,thresh)
comb <- comb[-4,]
colnames(comb) <- c('node_size','thresh')
n_cores <- 16

lapply(1:nrow(comb),function(i) RunBetaErr(resultPath,paste0('0.75_',comb$node_size[i],'_',as.character(comb$thresh[i]),'/'),as.numeric(as.character(comb$thresh[i])),n_cores,F))
lapply(1:nrow(comb),function(i) RunBetaErr(resultPath,paste0('0.75_',comb$node_size[i],'_',as.character(comb$thresh[i]),'/'),as.numeric(as.character(comb$thresh[i])),n_cores,T))


#Plots comparing training and testing set treatment effects
# data <- readRDS(paste0(resultPath,'data_p_','3e-05','.rds'))
# training_testing_set <- ExtractSubSample(data,readRDS('~/bsu_scratch/UKB_Data/training_set.rds'),readRDS('~/bsu_scratch/UKB_Data/test_set.rds'))
# testing_set <- training_testing_set$outofbag
# 
# test_tree_30000_result <- pbmclapply(1:5000,function(i){
#   test_tree <- readRDS(paste0(resultPath,'0.75_30000_3e-5/','tree',i,'.rds'))
#   CompareBeta(testing_set,test_tree$bootstrapPartyTree)
# },mc.cores = 16,ignore.interactive = T)
# saveRDS(test_tree_35000_result,file = 'bsu_scratch/Random_Forest/rs3821843_rs7340705_rs113210396_rs312487_rs11719824_rs3774530_rs3821856_sbp/test_tree_35000.rds')
# test_tree_30000_train <- unlist(lapply(test_tree_30000_result,function(x) x$trainingMainEffects))
# test_tree_30000_test <- unlist(lapply(test_tree_30000_result,function(x) x$testingMainEffects))
# # library(ggplot2)
# p1 <- ggplot(data.frame(training=test_tree_30000_train,testing=test_tree_30000_test), aes(x=training, y=testing)) + xlim(-1,1) + ylim(-1,1) +  geom_bin2d() + theme_bw() + labs(x='Training Set Treatment Effect',y='Testing Set Treatment Effect') + ggtitle('P<3e-5, MinNode=30000')
# # 
# # 
# 
# data <- readRDS(paste0(resultPath,'data_p_','1e-05','.rds'))
# training_testing_set <- ExtractSubSample(data,readRDS('~/bsu_scratch/UKB_Data/training_set.rds'),readRDS('~/bsu_scratch/UKB_Data/test_set.rds'))
# testing_set <- training_testing_set$outofbag
# 
# test_tree_10000_result <- pbmclapply(1:5000,function(i){
#   test_tree <- readRDS(paste0(resultPath,'0.75_10000_3e-5/','tree',i,'.rds'))
#   CompareBeta(testing_set,test_tree$bootstrapPartyTree)
# },mc.cores = 30,ignore.interactive = T)
# saveRDS(test_tree_10000_result,file = 'bsu_scratch/Random_Forest/rs3821843_rs7340705_rs113210396_rs312487_rs11719824_rs3774530_rs3821856_sbp/test_tree_10000.rds')
# test_tree_10000_train <- unlist(lapply(test_tree_10000_result,function(x) x$trainingMainEffects))
# test_tree_10000_test <- unlist(lapply(test_tree_10000_result,function(x) x$testingMainEffects))
# library(ggplot2)
# p2 <- ggplot(data.frame(training=test_tree_10000_train,testing=test_tree_10000_test), aes(x=training, y=testing)) + geom_bin2d() + theme_bw() + labs(x='Training Set Treatment Effect',y='Testing Set Treatment Effect') + ggtitle('P<3e-5, MinNode=10000')
# library(ggpubr)
# ggarrange(plotlist = list(p1,p2),ncol = 1,nrow = 2)
