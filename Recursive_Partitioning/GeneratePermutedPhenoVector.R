####################################################################################
#Generate a list of permuted phenotype vectors 
#Input: p_value treshold, stored unpermuted data, and number of permutations.
#Output: .rds file containing the list of genotype matrices.
####################################################################################

source('~/MRC_BSU_Internship/Recursive_Partitioning/ExtractSubsample.R')
#Permuted dosgage of each SNP, across samples
GeneratePermutation <- function(phenoVector){
  permPheno <- sample(phenoVector,size = length(phenoVector),replace = F)
  return(permPheno)
}

#Generate permuted genotype matrix
GeneratePermutatedPhenoMatrix <- function(resultPath,suffix,p_thresh,n_perm){
  #Read in unpermuted data, and extract testing set
  data <- readRDS(file = paste0(resultPath,'data_p_',p_thresh,'.rds'))
  training_testing_set <- ExtractSubSample(data,readRDS('~/bsu_scratch/UKB_Data/training_set.rds'),readRDS('~/bsu_scratch/UKB_Data/test_set.rds'))
  testingSetSamples <- training_testing_set$outofbag
  
  #Extract genotype matrix from testing set
  rownames(testingSetSamples$dosageMatrix) <- NULL
  names(testingSetSamples$dosageTarget) <- NULL
  names(testingSetSamples$phenotypes) <- NULL
  rownames(testingSetSamples$covariates) <- NULL
  
  #Generate permuted genotype matrix, 
  permPhenoVector <- lapply(1:n_perm,function(x) GeneratePermutation(testingSetSamples$phenotypes))
  saveRDS(permPhenoVector,file = paste0(paste0(resultPath,'perm_phenotype_p_',p_thresh,'.rds')))
}
#args=(commandArgs(TRUE))
# thresh <- args[[1]]
thresh <- c('1e-5','2e-5','3e-5','4e-5')
lapply(thresh,function(x) GeneratePermutatedPhenoMatrix(resultPath = '~/bsu_scratch/Random_Forest/rs3821843_rs7340705_rs113210396_rs312487_rs11719824_rs3774530_rs3821856_sbp/',
                             suffix = paste0('0.75_',node_size,'_',thresh,'/'),p_thresh = as.numeric(x),n_perm = 1000))
