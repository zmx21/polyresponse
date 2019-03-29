####################################################################################
#Generate a list of permuted genotype matrix (for use as predictor in RF). 
# permuted dosgage of each SNP indepently, across samples. This maintains MAF of each SNP.
#Input: p_value treshold, stored unpermuted data, and number of permutations.
#Output: .rds file containing the list of genotype matrices.
####################################################################################

source('~/MRC_BSU_Internship/Recursive_Partitioning/ExtractSubsample.R')
#Permuted dosgage of each SNP, across samples
GeneratePermutation <- function(genotypeMatrix){
  return(as.data.frame(apply(genotypeMatrix,2,function(x) sample(x,size = length(x),replace = F))))
}

#Generate permuted genotype matrix
GeneratePermutatedGenoMatrix <- function(resultPath,suffix,p_thresh,n_perm){
  #Read in unpermuted data, and extract testing set
  data <- readRDS(file = paste0(resultPath,'data_p_',p_thresh,'.rds'))
  training_testing_set <- ExtractSubSample(data,readRDS('~/bsu_scratch/UKB_Data/training_set.rds'),readRDS('~/bsu_scratch/UKB_Data/test_set.rds'))
  testingSetSamples <- training_testing_set$outofbag
  
  #Extract genotype matrix from testing set
  genotypeMatrix <- testingSetSamples$dosageMatrix
  rownames(genotypeMatrix) <- NULL
  names(testingSetSamples$dosageTarget) <- NULL
  names(testingSetSamples$phenotypes) <- NULL
  rownames(testingSetSamples$covariates) <- NULL
  
  #Generate permuted genotype matrix, 
  permGenotypeMatrix <- lapply(1:n_perm,function(x) GeneratePermutation(genotypeMatrix))
  saveRDS(permGenotypeMatrix,file = paste0(paste0(resultPath,'perm_genotype_p_',p_thresh,'.rds')))
}
args=(commandArgs(TRUE))
thresh <- args[[1]]
GeneratePermutatedGenoMatrix(resultPath = '~/bsu_scratch/Random_Forest/rs3821843_rs7340705_rs113210396_rs312487_rs11719824_rs3774530_rs3821856_sbp/',
                             suffix = paste0('0.75_',node_size,'_',thresh,'/'),p_thresh = as.numeric(thresh),n_perm = 1000)