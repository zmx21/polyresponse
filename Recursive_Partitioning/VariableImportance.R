####################################################################################
#Appens variable improtance measures to each tree in RF (already constructed)
#Input: Line 39
#Output: saved .rds file of trees, with the same name but with variable importance info appended.
####################################################################################


source('~/MRC_BSU_Internship/Recursive_Partitioning/InteractionTree.R')
source('~/MRC_BSU_Internship/Recursive_Partitioning/ExtractSubsample.R')
library(pbmcapply)
library(parallel)
library(partykit)
VariableImportance <- function(bootstrapTree,outOfBagData){
  #Send out of bag samples down tree to calculate total interaction
  totalInteraction <- sum(CalculateTotalInteraction(as.data.frame(outOfBagData$dosageMatrix),bootstrapTree,outOfBagData$dosageTarget,outOfBagData$phenotypes,outOfBagData$covariates))
  #loop through all features, calculate invidual variable importance
  allNodeInfo <- unlist(nodeapply(bootstrapTree,ids = nodeids(bootstrapTree),FUN=function(x) x$info))
  inclFeatures <- unique(allNodeInfo[sapply(allNodeInfo,function(x) grepl('rs',x))])
  featureInteractions <- rep(NA,ncol(outOfBagData$dosageMatrix))
  names(featureInteractions) <- colnames(outOfBagData$dosageMatrix)
  for(i in 1:length(featureInteractions)){
    if(!names(featureInteractions)[i]%in%inclFeatures){
      next
    }
    permutedDosageMatrix <- outOfBagData$dosageMatrix
    permutedDosageMatrix[,i] <- sample(as.vector(permutedDosageMatrix[,i]),size = length(as.vector(permutedDosageMatrix[,i])),replace = F)
    featureInteractions[i] <- sum(CalculateTotalInteraction(as.data.frame(permutedDosageMatrix),bootstrapTree,outOfBagData$dosageTarget,outOfBagData$phenotypes,outOfBagData$covariates))
  }
  variableImportance <- (totalInteraction - featureInteractions) / totalInteraction
  return(variableImportance)
}
args=(commandArgs(TRUE))
if(length(args)==0){
  resultPath <- '~/bsu_scratch/Random_Forest/rs3821843_rs7340705_rs113210396_rs312487_rs11719824_rs3774530_rs3821856_sbp/'
  p_thresh <-  '3e-5'
  node_size <- 35000
  n_cores <- 16
}else if (length(args) != 4){
  stop("You need to supply:\n",
       "# 1: RF result path\n",
       "# 2: P-value threshold\n",
       "# 3: Minimum node size\n",
       "# 4: Number of cores",
       "Exiting...", call.=FALSE)
}else{
  print('All Arguments Supplied')
  str <- c("# 1: RF result path:",
           '# 2: P-value threshold:',
           "# 3: Minimum node size:",
           "# 4: Number of cores:")
  variables <- c('resultPath','p_thresh','node_size','n_cores')
  for(i in 1:length(args)){
    eval(parse(text=paste0(variables[i],'=',"'",args[[i]],"'")))
    print(paste0(str[i],"   ",args[i]))
  }
}
RunVarImp <- function(treeIndex,resultPath,node_size,p_thresh,testing_set,training_set){
  currentTree <- readRDS(paste0(resultPath,'0.75_',node_size,'_',p_thresh,'/tree',treeIndex,'.rds'))
  currentTree$testingVarImp <- VariableImportance(currentTree$bootstrapPartyTree,testing_set)
  currentTree$trainingVarImp <- VariableImportance(currentTree$bootstrapPartyTree,ExtractSubSample(training_set,currentTree$bootstrapIndex,currentTree$outofbagIndex)$bootstrap)
  if(!is.null(currentTree$variableImportance)){
    currentTree$outOfBagVarImp <- currentTree$variableImportance
    currentTree$variableImportance <- NULL
  }
  saveRDS(currentTree,paste0(resultPath,'0.75_',node_size,'_',p_thresh,'/tree',treeIndex,'.rds'))
}
data <- readRDS(paste0(resultPath,'data_p_',as.numeric(p_thresh),'.rds'))
training_testing_set <- ExtractSubSample(data,readRDS('~/bsu_scratch/UKB_Data/training_set.rds'),readRDS('~/bsu_scratch/UKB_Data/test_set.rds'))
testing_set <- training_testing_set$outofbag
training_set <- training_testing_set$bootstrap
trash <- pbmclapply(1:5000,function(i) RunVarImp(i,resultPath,node_size,p_thresh,testing_set,training_set),ignore.interactive = T,mc.cores = n_cores)