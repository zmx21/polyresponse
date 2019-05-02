####################################################################################
#Calculates training set treatment effects, and calcualtes the Root Squared error (against training set)
#Input: path to where random forests and dosage data is stored,the thresh used to construct the RF, and 
#whether the true or permuted testing set should be used. 
#Output: list of numerics (root squared error)
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
  #Calculate difference between training and testing (error)
  err <- testingMainEffects-trainingMainEffects
  names(err) <- uniqueLeafNodes
  return(list(err=err,trainingMainEffects=trainingMainEffects,testingMainEffects=testingMainEffects))
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
  err <- lapply(perm_results,function(x) x-trainingMainEffects)
  return(err)
}
#Get treatment effects for all trees within random forest and collapse.
CalculateBetaError <- function(resultPath,suffix,p_thresh,n_cores,perm,chunks){
  chunks <- as.numeric(strsplit(chunks,':')[[1]][1]):as.numeric(strsplit(chunks,':')[[1]][2])
  if(perm=='Geno'){
    err <- pbmclapply(chunks,function(i){
      #Perm ordering is different
      #find all individual trees in random forest
      treePaths <- dir(paste0(resultPath,suffix))
      treePaths <- treePaths[sapply(treePaths,function(x) grepl('tree',x))]
      
      #Precalculated list of permuted predictions should be saved (saved by PermutedPrediction.R)
      perm_result <- readRDS(paste0(resultPath,suffix,'prediction_betas_perm/testing_main_effects_tree',which(treePaths==paste0('tree',i,'.rds')),'.rds'))
      test_tree <- readRDS(paste0(resultPath,suffix,'tree',i,'.rds'))$bootstrapPartyTree
      CompareBetaPerm(perm_result,test_tree)
    },mc.cores = n_cores,ignore.interactive = T)
    return(err)
  }else if(perm == 'None'){
    data <- readRDS(paste0(resultPath,'data_p_',p_thresh,'.rds'))
    training_testing_set <- ExtractSubSample(data,readRDS('~/bsu_scratch/LDL_Project_Data/Genotype_Data/training_set.rds'),
                                             readRDS('~/bsu_scratch/LDL_Project_Data/Genotype_Data/test_set.rds'))
    testing_set <- training_testing_set$outofbag
    err <- pbmclapply(chunks,function(i){
      test_tree <- readRDS(paste0(resultPath,suffix,'tree',i,'.rds'))
      CompareBeta(testing_set,test_tree$bootstrapPartyTree)$err
    },mc.cores = n_cores,ignore.interactive = T)
    return(err)
  }else{
    err <- pbmclapply(chunks,function(i){
      #Precalculated list of permuted predictions should be saved (saved by PermutedPrediction.R)
      perm_result <- readRDS(paste0(resultPath,suffix,'prediction_betas_pheno_perm/testing_main_effects_tree',i,'.rds'))
      test_tree <- readRDS(paste0(resultPath,suffix,'tree',i,'.rds'))$bootstrapPartyTree
      CompareBetaPerm(perm_result,test_tree)
    },mc.cores = n_cores,ignore.interactive = T)
    return(err)
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
    saveRDS(betaErr,paste0(resultPath,suffix,'beta_err',file_suffix,'.rds'))
  }else if(perm == 'Geno'){
    betaErr <- CalculateBetaError(resultPath,suffix,p_thresh,n_cores,perm,chunks)
    saveRDS(betaErr,paste0(resultPath,suffix,'beta_err_perm',file_suffix,'.rds'))
  }else{
    betaErr <- CalculateBetaError(resultPath,suffix,p_thresh,n_cores,perm,chunks)
    saveRDS(betaErr,paste0(resultPath,suffix,'beta_err_pheno_perm',file_suffix,'.rds'))
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
node_size=c(1000,2500,5000,10000,20000,30000,40000)
thresh <- c('5e-6','1e-5','3e-5','5e-5','7e-5','1e-4')
comb <- expand.grid(node_size,thresh)
colnames(comb) <- c('node_size','thresh')
n_cores <- 12
for(i in 1:nrow(comb)){
  print(comb[i,])
  RunBetaErr(resultPath,paste0('0.75_',comb$node_size[i],'_',as.character(comb$thresh[i]),'/'),as.numeric(as.character(comb$thresh[i])),n_cores,'Pheno',"1:500")
  RunBetaErr(resultPath,paste0('0.75_',comb$node_size[i],'_',as.character(comb$thresh[i]),'/'),as.numeric(as.character(comb$thresh[i])),n_cores,'Pheno',"501:1000")
  RunBetaErr(resultPath,paste0('0.75_',comb$node_size[i],'_',as.character(comb$thresh[i]),'/'),as.numeric(as.character(comb$thresh[i])),n_cores,'Pheno',"1001:1500")
  RunBetaErr(resultPath,paste0('0.75_',comb$node_size[i],'_',as.character(comb$thresh[i]),'/'),as.numeric(as.character(comb$thresh[i])),n_cores,'Pheno',"1501:2000")
  RunBetaErr(resultPath,paste0('0.75_',comb$node_size[i],'_',as.character(comb$thresh[i]),'/'),as.numeric(as.character(comb$thresh[i])),n_cores,'None',"1:500")
  RunBetaErr(resultPath,paste0('0.75_',comb$node_size[i],'_',as.character(comb$thresh[i]),'/'),as.numeric(as.character(comb$thresh[i])),n_cores,'None',"501:1000")
  RunBetaErr(resultPath,paste0('0.75_',comb$node_size[i],'_',as.character(comb$thresh[i]),'/'),as.numeric(as.character(comb$thresh[i])),n_cores,'None',"1001:1500")
  RunBetaErr(resultPath,paste0('0.75_',comb$node_size[i],'_',as.character(comb$thresh[i]),'/'),as.numeric(as.character(comb$thresh[i])),n_cores,'None',"1501:2000")
}
