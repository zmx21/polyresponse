source('~/MRC_BSU_Internship/Recursive_Partitioning/InteractionTree.R')
source('~/MRC_BSU_Internship/Recursive_Partitioning/ExtractSubsample.R')
source('~/MRC_BSU_Internship/Recursive_Partitioning/LoadDosage.R')
library(partykit)
library(pbmcapply)
CompareBeta <- function(testing_set,tree,perm){
  #Assign each sample within the testing set to a leaf node
  leafNodeAssignment <- predict(tree,as.data.frame(testing_set$dosageMatrix))
  uniqueLeafNodes <- unique(leafNodeAssignment)
  
  #If permutation is specified, permute the leaf node assignment
  if(perm){
    leafNodeAssignmentPerm <- leafNodeAssignment[sample(1:length(leafNodeAssignment),size = length(leafNodeAssignment),replace = F)]
    names(leafNodeAssignmentPerm) <- names(leafNodeAssignment)
    leafNodeAssignment <- leafNodeAssignmentPerm
  }
  #Calculate beta coefficients of the subgroups within the testing set.
  testingMainEffects <- rep(NA,length(uniqueLeafNodes))
  names(testingMainEffects) <- uniqueLeafNodes
  for(i in 1:length(uniqueLeafNodes)){
    currentIndicator <- leafNodeAssignment == uniqueLeafNodes[i]
    testingMainEffects[i] <- FitMainEffectModel(testing_set$dosageTarget[currentIndicator],testing_set$phenotypes[currentIndicator],testing_set$covariates[currentIndicator])$coeff['treatment']
  }
  trainingMainEffects <- lapply(uniqueLeafNodes,function(x) nodeapply(tree,ids = x,FUN = function(n) n$info))
  trainingMainEffects <- sapply(trainingMainEffects,function(x) as.numeric(unlist(strsplit(x = unlist(x),split = "  "))[1]))
  names(trainingMainEffects) <- uniqueLeafNodes
  rse <- sqrt((testingMainEffects-trainingMainEffects)^2)
  names(rse) <- NULL
  return(list(testingMainEffects=testingMainEffects,trainingMainEffects=trainingMainEffects,rse=rse))
}
CalculateBetaError <- function(resultPath,suffix,p_thresh,n_cores,perm){
  data <- readRDS(paste0(resultPath,'data_p_',p_thresh,'.rds'))
  training_testing_set <- ExtractSubSample(data,readRDS('~/bsu_scratch/UKB_Data/training_set.rds'),readRDS('~/bsu_scratch/UKB_Data/test_set.rds'))
  testing_set <- training_testing_set$outofbag
  rse <- pbmclapply(1:5000,function(i){
    test_tree <- readRDS(paste0(resultPath,suffix,'tree',i,'.rds'))
    CompareBeta(testing_set,test_tree$bootstrapPartyTree,perm)$rse
  },mc.cores = n_cores,ignore.interactive = T)
  return(do.call(c,rse))
}

RunBetaErr <- function(resultPath,suffix,p_thresh,n_cores,perm){
  for(i in 1:length(suffix)){
    if(!perm){
      betaErr <- CalculateBetaError(resultPath,suffix[i],p_thresh,n_cores,perm)
      saveRDS(betaErr,paste0(resultPath,suffix[i],'beta_err.rds'))
    }else{
      betaErr <- CalculateBetaError(resultPath,suffix[i],p_thresh,n_cores,perm)
      saveRDS(betaErr,paste0(resultPath,suffix[i],'beta_err_perm.rds'))
    }
  }
}
args=(commandArgs(TRUE))
node_size <- as.numeric(args[[1]])
n_cores <- as.numeric(args[[2]])
resultPath <- '~/bsu_scratch/Random_Forest_0p75/rs3821843_rs7340705_rs113210396_rs312487_rs11719824_rs3774530_rs3821856_sbp/'
thresh <- c('1e-5','2e-5','3e-5','4e-5','6e-5')
for(i in 1:length(thresh)){
  RunBetaErr(resultPath,paste0('0.75_',node_size,'_',thresh[i],'/'),as.numeric(thresh[i]),n_cores,T)
  RunBetaErr(resultPath,paste0('0.75_',node_size,'_',thresh[i],'/'),as.numeric(thresh[i]),n_cores,F)
}

# data <- readRDS(paste0(resultPath,'data_p_','3e-05','.rds'))
# training_testing_set <- ExtractSubSample(data,readRDS('~/bsu_scratch/UKB_Data/training_set.rds'),readRDS('~/bsu_scratch/UKB_Data/test_set.rds'))
# testing_set <- training_testing_set$outofbag

# test_tree_35000_result <- pbmclapply(1:5000,function(i){
#   test_tree <- readRDS(paste0(resultPath,'0.75_35000_3e-5/','tree',i,'.rds'))
#   CompareBeta(testing_set,test_tree$bootstrapPartyTree,F)
# },mc.cores = 30,ignore.interactive = T)
# saveRDS(test_tree_35000_result,file = 'bsu_scratch/Random_Forest/rs3821843_rs7340705_rs113210396_rs312487_rs11719824_rs3774530_rs3821856_sbp/test_tree_35000.rds')
# test_tree_35000_train <- unlist(lapply(test_tree_35000_result,function(x) x$trainingMainEffects))
# test_tree_35000_test <- unlist(lapply(test_tree_35000_result,function(x) x$testingMainEffects))
# library(ggplot2)
# p1 <- ggplot(data.frame(training=test_tree_35000_train,testing=test_tree_35000_test), aes(x=training, y=testing)) + xlim(-1,1) + ylim(-1,1) +  geom_bin2d() + theme_bw() + labs(x='Training Set Treatment Effect',y='Testing Set Treatment Effect') + ggtitle('P<3e-5, MinNode=35000')
# 
# 
# test_tree_10000_result <- pbmclapply(1:5000,function(i){
#   test_tree <- readRDS(paste0(resultPath,'0.75_10000_3e-5/','tree',i,'.rds'))
#   CompareBeta(testing_set,test_tree$bootstrapPartyTree,F)
# },mc.cores = 30,ignore.interactive = T)
# saveRDS(test_tree_10000_result,file = 'bsu_scratch/Random_Forest/rs3821843_rs7340705_rs113210396_rs312487_rs11719824_rs3774530_rs3821856_sbp/test_tree_10000.rds')
# test_tree_10000_train <- unlist(lapply(test_tree_10000_result,function(x) x$trainingMainEffects))
# test_tree_10000_test <- unlist(lapply(test_tree_10000_result,function(x) x$testingMainEffects))
# library(ggplot2)
# p2 <- ggplot(data.frame(training=test_tree_10000_train,testing=test_tree_10000_test), aes(x=training, y=testing)) + geom_bin2d() + theme_bw() + labs(x='Training Set Treatment Effect',y='Testing Set Treatment Effect') + ggtitle('P<3e-5, MinNode=10000')
# library(ggpubr)
# ggarrange(plotlist = list(p1,p2),ncol = 1,nrow = 2)
