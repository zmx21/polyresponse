####################################################################################
#Seperate UKB data into a training and testing set. Stores indices as RDS files.
####################################################################################
source('../Load_Phenotype/Load_Phenotype.R')
GenerateTrainingTestSet <- function(splitFraction){
  path <-  '~/bsu_scratch/LDL_Project_Data/Genotype_Data/'
  sample_file_prefix <- 'ukbb_LDL_metadata_with_PC'
  phenotype = 'LDLdirect'
  cov_names <- c('sex','ages','bmi')
  eur_only = 1
  med = 1
  
  #Load Phenotype and Covariates
  phenotypesAndCov <- LoadPhenotype(path,sample_file_prefix,phenotype,cov_names,eur_only,med)
  
  #Training set 
  training_samples <- sample(1:sum(phenotypesAndCov$samplesToKeep),size = floor(splitFraction * sum(phenotypesAndCov$samplesToKeep)),replace = F)
  test_samples <- setdiff(1:sum(phenotypesAndCov$samplesToKeep),training_samples)
  saveRDS(training_samples,paste0(path,'training_set.rds'))
  saveRDS(test_samples,paste0(path,'test_set.rds')) 
}
