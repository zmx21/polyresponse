source('~/MRC_BSU_Internship_LDL/Recursive_Partitioning/ExtractSubsample.R')
source('~/MRC_BSU_Internship_LDL/Recursive_Partitioning/InteractionTree.R')

library(partykit)
library(parallel)
library(pbmcapply)
#Runs predictions of permuted samples in chunks
PermutedPhenoInteraction <- function(testingSetSamples,randomForestPath,n_cores,tree_chunks,permPhenotypeVect){
  system(paste0('mkdir -p ',randomForestPath,'/interactions_pheno_perm/'))
  #Parse tree chunk argument
  tree_chunks <- unlist(strsplit(tree_chunks,':'))
  tree_chunks <- tree_chunks[1]:tree_chunks[2]
  
  res <- pbmclapply(tree_chunks,function(i) {
  #res <- lapply(tree_chunks,function(i){
    curTree <- readRDS(paste0(randomForestPath,'/tree',i,'.rds'))
    #Get Total interaction of testing set based tree, using random phenotype
    permPhenotypeVect <- lapply(permPhenotypeVect,function(x){
      names(x) <- names(testingSetSamples$phenotypes);
      return(x)
    })
    interactions <- pbmclapply(permPhenotypeVect[1:50],function(x) CalculateTotalInteraction(as.data.frame(testingSetSamples$dosageMatrix),
                                                                                   curTree$bootstrapPartyTree,testingSetSamples$dosageTarget,
                                                                                   x,testingSetSamples$covariates),mc.cores = n_cores,ignore.interactive = F)
    saveRDS(interactions,paste0(randomForestPath,'/interactions_pheno_perm/','interactions_tree',i,'.rds'))
  },mc.cores = 1,ignore.interactive = T)
  #})
}
#MAIN FUNCTION
RunPermutedPhenoInteraction <- function(resultPath,suffix,p_thresh,n_cores,tree_chunks){
  #Load training set samples
  data <- readRDS(paste0(resultPath,'data_p_',p_thresh,'.rds'))
  training_testing_set <- ExtractSubSample(data,readRDS('~/bsu_scratch/LDL_Project_Data/Genotype_Data/training_set.rds'),
                                           readRDS('~/bsu_scratch/LDL_Project_Data/Genotype_Data/test_set.rds'))
  remove(data);gc(verbose = F);
  testing_set <- training_testing_set$outofbag
  remove(training_testing_set);gc(verbose = F);
  #Load permuted phenotype vectors
  permPhenotypeVect <- readRDS(paste0(resultPath,'perm_phenotype_p_',p_thresh,'.rds'))
  PermutedPhenoInteraction(testing_set,paste0(resultPath,suffix),n_cores,tree_chunks,permPhenotypeVect)
}
args=(commandArgs(TRUE))
node_size <- as.numeric(args[[1]])
thresh <- args[[2]]
n_cores <- as.numeric(args[[3]])
tree_chunks <- args[[4]]
# node_size <- 40000
# thresh <- '5e-6'
# n_cores <- 1
# tree_chunks <- '1:1'
print(c('node_size'=node_size,'thresh'=thresh,'n_cores'=n_cores,'tree_chunks'=tree_chunks))
resultPath <- '~/bsu_scratch/LDL_Project_Data/Random_Forest/rs12916_rs17238484_rs5909_rs2303152_rs10066707_rs2006760_LDLdirect/'
#resultPath <- '~/bsu_scratch/LDL_Project_Data/Random_Forest/rs12916_rs72633963_rs55727654_rs17648121_rs2303152_rs62366588_rs75240579_rs111353455_LDLdirect/'

RunPermutedPhenoInteraction(resultPath,paste0('0.75_',node_size,'_',thresh,'/'),as.numeric(thresh),n_cores,tree_chunks)