####################################################################################
#Calculates training set treatment effects, and calcualtes the weighted standard deviation
#Input: path to where random forests and dosage data is stored,the thresh used to construct the RF, and 
#whether the true or permuted testing set should be used. 
#Output: list of numerics (weighted sd)
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
  #randAssignment <- lapply(1:1000,function(x) sample(leafNodeAssignment,size = length(leafNodeAssignment),replace = F))
  uniqueLeafNodes <- sort(unique(leafNodeAssignment),decreasing = F)
  
  #Calculate beta coefficients of the subgroups within the testing set.
  testingBeta <- rep(NA,length(uniqueLeafNodes))
  #randBeta <- matrix(NaN,nrow = length(uniqueLeafNodes),ncol=1000)
  
  testingBetaFit <- list()
  #randBetaFit <- list()
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
  
  #Calculate weighted sd of training and testing set
  sdTraining <- sqrt(sum(trainingSampleSizes * (trainingBeta - weighted.mean(trainingBeta,trainingSampleSizes))^2)/sum(trainingSampleSizes))
  sdTesting <- sqrt(sum(testingSampleSizes * (testingBeta - weighted.mean(testingBeta,testingSampleSizes))^2)/sum(testingSampleSizes))
  return(list(sdTraining=sdTraining,sdTesting=sdTesting))
}
#Compare training set and PERMUTED testing set treatment effects. 
#Precalculated list of permuted predictions should be provided (saved by PermutedPrediction.R)
CompareBetaPerm <- function(perm_results,perm_node_size){
  sdPerm <- lapply(perm_results,function(x) sqrt(sum(perm_node_size * (x - weighted.mean(x,perm_node_size))^2)/sum(perm_node_size)))
  return(sdPerm)
}
#Get treatment effects for all trees within random forest and collapse.
CalculateBetaError <- function(resultPath,suffix,p_thresh,n_cores,perm,chunks,MAF){
  chunks <- as.numeric(strsplit(chunks,':')[[1]][1]):as.numeric(strsplit(chunks,':')[[1]][2])
  if(perm == 'None'){
    data <- readRDS(paste0(resultPath,'data_p_',p_thresh,'_maf_',MAF,'.rds'))
    training_testing_set <- ExtractSubSample(data,readRDS('~/bsu_scratch/LDL_Project_Data_Aug2019/Genotype_Data/training_set.rds'),
                                             readRDS('~/bsu_scratch/LDL_Project_Data_Aug2019/Genotype_Data/test_set.rds'))
    testing_set <- training_testing_set$outofbag
    sd <- pbmclapply(chunks,function(i){
      test_tree <- readRDS(paste0(resultPath,suffix,'tree',i,'.rds'))
      CompareBeta(testing_set,test_tree$bootstrapPartyTree)
    #})
    },mc.cores = n_cores,ignore.interactive = T)
    return(sd)
  }
  else if(perm=='Pheno'){
    sd <- pbmclapply(chunks,function(i){ #lapply(chunks,function(i){
      #Precalculated list of permuted predictions should be saved (saved by PermutedPrediction.R)
      perm_result <- readRDS(paste0(resultPath,suffix,'prediction_betas_perm/testing_main_effects_tree',i,'.rds'))
      perm_assignment <-readRDS(paste0(resultPath,suffix,'prediction_betas_perm/testing_indiv_main_effects_tree',i,'.rds'))
      #Find sample size in each node
      perm_node_size <- table(names(perm_result[[1]])[match(perm_assignment[[1]],perm_result[[1]])])
      #reorder the nodes
      perm_node_size <- perm_node_size[match(names(perm_result[[1]]),names(perm_node_size))]
      return(CompareBetaPerm(perm_result,perm_node_size))
    #})
    },mc.cores = n_cores,ignore.interactive = T)
    return(sd)
  }
}

#MAIN FUNCTION. If non-permuted testing set is to be used, calculation is on the fly. 
#If permuted set to be used, precalculated list of permuted predictions should be saved (saved by PermutedPrediction.R)

RunBetaErr <- function(resultPath,suffix,p_thresh,n_cores,perm,chunks='1:2000',MAF){
  file_suffix <- ''
  if(chunks != '1:2000'){
    file_suffix <- paste0('_',strsplit(chunks,':')[[1]][1],'_',strsplit(chunks,':')[[1]][2])
  }
  if(perm=='None'){
    betaErr <- CalculateBetaError(resultPath,suffix,p_thresh,n_cores,perm,chunks,MAF)
    saveRDS(betaErr,paste0(resultPath,suffix,'sd',file_suffix,'.rds'))
  }else if (perm =='Pheno'){
    betaErr <- CalculateBetaError(resultPath,suffix,p_thresh,n_cores,perm,chunks,MAF)
    saveRDS(betaErr,paste0(resultPath,suffix,'sd_pheno_perm_new',file_suffix,'.rds'))
  }
}
resultPath <- '~/bsu_scratch/LDL_Project_Data_Aug2019/Random_Forest_Old/rs12916_rs17238484_rs5909_rs2303152_rs10066707_rs2006760_LDLdirect/'
# args=(commandArgs(TRUE))
# node_size <- as.numeric(args[[1]])
# thresh <- args[[2]]
# n_cores <- as.numeric(args[[3]])
# MAF = args[[4]]
# node_size= c(10000,20000,30000,40000) #c(1000,2500,5000,10000,20000,30000,40000)
# thresh <- c('6.5e-6','7e-6','7.5e-6','8e-6','9e-5') #c('5e-6','1e-5','3e-5','5e-5','7e-5','1e-4')
MAF = '5e-2'
# comb <- expand.grid('node_size'=node_size,'thresh'=thresh)

node_size <- c(10000,20000,30000,40000)
thresh <- c('6.75e-6','6.5e-6','6e-6','5e-6','3e-6') #c('9e-6','7.5e-6','7.25e-6','7e-6','6.75e-6','6.5e-6','6e-6','5e-6','3e-6')
comb <- expand.grid('node_size'=node_size,'thresh'=thresh)
thresh <- c('9e-5','7e-5','5e-5','3e-5','1e-5')
node_size <- c(5000,10000,20000,30000,40000)
comb <- rbind(comb,expand.grid('node_size'=node_size,'thresh'=thresh))


n_cores <- 20
for(i in 1:nrow(comb)){
  print(comb[i,])
  # RunBetaErr(resultPath,paste0('0.75_',comb$node_size[i],'_',as.character(comb$thresh[i]),'_',MAF,'/'),as.numeric(as.character(comb$thresh[i])),n_cores,'None',"1:500",MAF)
  # RunBetaErr(resultPath,paste0('0.75_',comb$node_size[i],'_',as.character(comb$thresh[i]),'_',MAF,'/'),as.numeric(as.character(comb$thresh[i])),n_cores,'None',"501:1000",MAF)
  # RunBetaErr(resultPath,paste0('0.75_',comb$node_size[i],'_',as.character(comb$thresh[i]),'_',MAF,'/'),as.numeric(as.character(comb$thresh[i])),n_cores,'None',"1001:1500",MAF)
  # RunBetaErr(resultPath,paste0('0.75_',comb$node_size[i],'_',as.character(comb$thresh[i]),'_',MAF,'/'),as.numeric(as.character(comb$thresh[i])),n_cores,'None',"1501:2000",MAF)
  RunBetaErr(resultPath,paste0('0.75_',comb$node_size[i],'_',as.character(comb$thresh[i]),'_',MAF,'/'),as.numeric(as.character(comb$thresh[i])),n_cores,'Pheno',"1:500",MAF)
  RunBetaErr(resultPath,paste0('0.75_',comb$node_size[i],'_',as.character(comb$thresh[i]),'_',MAF,'/'),as.numeric(as.character(comb$thresh[i])),n_cores,'Pheno',"501:1000",MAF)
  RunBetaErr(resultPath,paste0('0.75_',comb$node_size[i],'_',as.character(comb$thresh[i]),'_',MAF,'/'),as.numeric(as.character(comb$thresh[i])),n_cores,'Pheno',"1001:1500",MAF)
  RunBetaErr(resultPath,paste0('0.75_',comb$node_size[i],'_',as.character(comb$thresh[i]),'_',MAF,'/'),as.numeric(as.character(comb$thresh[i])),n_cores,'Pheno',"1501:2000",MAF)
}
