####################################################################################
#Generate a list of permuted phenotype vectors 
#Input: p_value treshold, stored unpermuted data, and number of permutations.
#Output: .rds file containing the list of genotype matrices.
####################################################################################

source('~/MRC_BSU_Internship_LDL/Recursive_Partitioning/ExtractSubsample.R')
#Permuted dosgage of each SNP, across samples
GeneratePermutation <- function(phenoVector){
  permPheno <- sample(phenoVector,size = length(phenoVector),replace = F)
  return(permPheno)
}

#Generate permuted genotype matrix
GeneratePermutatedPhenoMatrix <- function(resultPath,p_thresh,n_perm){
  #Read in unpermuted data, and extract testing set
  print(paste0(resultPath,'data_p_',p_thresh,'.rds'))
  data <- readRDS(file = paste0(resultPath,'data_p_',p_thresh,'.rds'))
  training_testing_set <- ExtractSubSample(data,
                                           readRDS('~/bsu_scratch/LDL_Project_Data/Genotype_Data/training_set.rds'),
                                           readRDS('~/bsu_scratch/LDL_Project_Data/Genotype_Data/test_set.rds'))
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
#thresh <- args[[1]]
#thresh <- c('5e-6','1e-5','3e-5','5e-6')
thresh <- '7e-5'
resultPath <- '~/bsu_scratch/LDL_Project_Data/Random_Forest/rs12916_rs17238484_rs5909_rs2303152_rs10066707_rs2006760_LDLdirect/'
#resultPath <- '~/bsu_scratch/LDL_Project_Data/Random_Forest/rs12916_rs72633963_rs55727654_rs17648121_rs2303152_rs62366588_rs75240579_rs111353455_LDLdirect/'

lapply(thresh,function(x) GeneratePermutatedPhenoMatrix(resultPath,p_thresh = as.numeric(x),n_perm = 1000))
