####################################################################################
#Apply RF to a list of permuted phenotype vector 
#Input: p_value treshold (of predictors to include), and minimum node size of tree.
#Output: .rds files containing predicted treatment effect of each leaf node
####################################################################################
source('~/MRC_BSU_Internship_LDL/Recursive_Partitioning/ExtractSubsample.R')
source('~/MRC_BSU_Internship_LDL/Recursive_Partitioning/InteractionTree.R')

library(partykit)
library(parallel)
library(pbmcapply)
#Generate training and testing set prediction, given a phenotype matrix
GenerateTestingPrediction <- function(tree,curPerm,testingSetData){
  nodeAssignment <- predict(tree,as.data.frame(testingSetData$dosageMatrix))
  nodeAssignment <- as.vector(nodeAssignment)
  uniqueLeafNodes <- sort(unique(nodeAssignment),decreasing = F)
  #Calculate treatment effect estimates based on testing set
  testingMainEffects <- rep(NA,length(uniqueLeafNodes))
  names(testingMainEffects) <- uniqueLeafNodes
  for(i in 1:length(uniqueLeafNodes)){
    currentIndicator <- nodeAssignment == uniqueLeafNodes[i]
    testingMainEffects[i] <- FitMainEffectModel(testingSetData$dosageTarget[curPerm][currentIndicator],
                                                testingSetData$phenotypes[curPerm][currentIndicator],
                                                testingSetData$covariates[curPerm,][currentIndicator,])$coeff['treatment']
  }
  mainEffectPred <- testingMainEffects[as.character(nodeAssignment)]
  names(mainEffectPred) <- NULL
  return(list(mainEffectPred=mainEffectPred,testingMainEffects=testingMainEffects))
}

#Runs predictions of permuted samples in chunks
PermutedPhenoPrediction <- function(testingSetData,randomForestPath,n_cores,tree_chunks,permVector){
  system(paste0('mkdir -p ',randomForestPath,'/prediction_betas_perm/'))
  
  #Parse tree chunk argument
  tree_chunks <- unlist(strsplit(tree_chunks,':'))
  tree_chunks <- tree_chunks[1]:tree_chunks[2]
  
  #Initialize summation matrix
  #Loop through all trees
  pb <- progressBar()
  for(i in 1:(length(tree_chunks))){
    setTxtProgressBar(pb,i/length(tree_chunks))
    #Load trees
    treeObj <- readRDS(paste0(randomForestPath,'tree',tree_chunks[i],'.rds'))
    #Predict treatment effects across bootstrap samples
    permResults <- lapply(permVector,function(k){
      GenerateTestingPrediction(treeObj$bootstrapPartyTree,
                                k,
                                testingSetData)
    })#,mc.cores = n_cores)
    testingMainEffects <- lapply(permResults,function(x) x$testingMainEffects)
    mainEffectPred <- lapply(permResults,function(x) x$mainEffectPred)
    saveRDS(testingMainEffects,file=paste0(randomForestPath,'/prediction_betas_perm/','testing_main_effects_tree',tree_chunks[i],'.rds'))
    saveRDS(mainEffectPred,file=paste0(randomForestPath,'/prediction_betas_perm/','testing_indiv_main_effects_tree',tree_chunks[i],'.rds'))
    }
}

#MAIN FUNCTION
RunPermutedPhenoPrediction <- function(resultPath,suffix,p_thresh,n_cores,tree_chunks,MAF){
  #Load training set samples
  data <- readRDS(paste0(resultPath,'data_p_',as.numeric(p_thresh),'_maf_',MAF,'.rds'))
  trainingTestingSet <- ExtractSubSample(data,readRDS('~/bsu_scratch/LDL_Project_Data/Genotype_Data/training_set.rds'),
                                           readRDS('~/bsu_scratch/LDL_Project_Data/Genotype_Data/test_set.rds'))
  remove(data);gc(verbose = F);
  testingSetData <- trainingTestingSet$outofbag
  remove(trainingTestingSet);gc(verbose = F);
  
  permVector <- readRDS(paste0(resultPath,'/randPerm.rds'))
  #Load permuted phenotype vectors
  PermutedPhenoPrediction(testingSetData,paste0(resultPath,suffix),n_cores,tree_chunks,permVector)
}

# args=(commandArgs(TRUE))
# node_size <- as.numeric(args[[1]])
# thresh <- args[[2]]
# n_cores <- as.numeric(args[[3]])
# tree_chunks <- args[[4]]
# MAF <- args[[5]]
# node_size <- 30000
# thresh <- '7e-6'
# n_cores <- 1
tree_chunks <- '1:2000'
MAF <- '5e-2'
# print(c('node_size'=node_size,'thresh'=thresh,'n_cores'=n_cores,'tree_chunks'=tree_chunks,'MAF'=MAF))
resultPath <- '~/bsu_scratch/LDL_Project_Data_Aug2019/Random_Forest_Old/12916_rs17238484_rs5909_rs2303152_rs10066707_rs2006760_LDLdirect/'
#resultPath <- '~/bsu_scratch/LDL_Project_Data/Random_Forest/rs12916_rs72633963_rs55727654_rs17648121_rs2303152_rs62366588_rs75240579_rs111353455_LDLdirect/'
#resultPath <- '~/bsu_scratch/LDL_Project_Data/Random_Forest/rs11206510_rs2479409_rs2149041_rs2479394_rs10888897_rs7552841_rs562556_LDLdirect/'
n_cores <- 10

node_size <- c(10000,20000,30000,40000)
thresh <- c('7.5e-6','7.25e-6','6.75e-6','6.5e-6','6e-6')
comb <- expand.grid(node_size,thresh)
node_size <- c(5000,node_size)
thresh <- c('9e-6','7e-6','5e-6','3e-6','9e-5','7e-5','5e-5','3e-5','1e-5')
comb <- rbind(comb,expand.grid(node_size,thresh))
colnames(comb) <- c('node_size','thresh')

print(c("nodesize"=node_size,"thresh"=thresh,"n_cores"=n_cores))
resultPath <- '~/bsu_scratch/LDL_Project_Data_Aug2019/Random_Forest_Old/rs12916_rs17238484_rs5909_rs2303152_rs10066707_rs2006760_LDLdirect/'

for(i in 1:nrow(comb)){
  RunPermutedPhenoPrediction(resultPath,paste0('0.75_',comb$node_size[i],'_',comb$thresh[i],'_5e-2/'),as.numeric(as.character(comb$thresh[i])),n_cores,tree_chunks,MAF)
}
