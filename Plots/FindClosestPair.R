library(dplyr)
library(data.table)
source('~/MRC_BSU_Internship_LDL/Recursive_Partitioning/ExtractSubsample.R')
cov_geno_data <- readRDS('~/bsu_scratch/LDL_Project_Data/Random_Forest/rs12916_rs17238484_rs5909_rs2303152_rs10066707_rs2006760_LDLdirect/data_p_1e-04.rds')
# int_db <- data.table::fread('~/bsu_scratch/LDL_Project_Data/Interaction_Data/HMGCR_LDL_known.txt')

training_testing_set <- ExtractSubSample(cov_geno_data,
                                         readRDS('~/bsu_scratch/LDL_Project_Data/Genotype_Data/training_set.rds'),
                                         readRDS('~/bsu_scratch/LDL_Project_Data/Genotype_Data/test_set.rds'))
testing_set <- training_testing_set$outofbag
training_set <- training_testing_set$bootstrap


TransformCovariantes <- function(cov){
  no_transform <- c()
  to_transfrom <- which(!colnames(cov) %in% no_transform)
  for(j in to_transfrom){
    cov[,j] <- (cov[,j] - mean(cov[,j])) / sd(cov[,j])
  }
  return(cov)
}
testing_set$covariates <- TransformCovariantes(testing_set$covariates)
training_set$covariates <- TransformCovariantes(training_set$covariates)

overall_fit <- lm(formula = paste0('LDL ~',paste(colnames(testing_set$covariates),collapse = '+')),data = cbind(data.frame(LDL = c(testing_set$phenotypes,training_set$phenotypes)),rbind(testing_set$covariates,training_set$covariates)))
fit_coeff <- overall_fit$coefficients[names(overall_fit$coefficients) != '(Intercept)']

FindCovDist <- function(training_set,testing_set,index,fit_coeff){
  training_cov <- training_set$covariates[index,]
  
  cov_dist <- apply(testing_set$covariates,1,function(x) as.numeric(as.matrix(training_cov - x) %*% as.matrix(fit_coeff)))
  return(cov_dist)
}

int_coeff <- int_db %>% dplyr::filter(rsid %in% colnames(training_set$dosageMatrix))
int_coeff <- -log10(int_coeff[match(int_coeff$rsid,colnames(training_set$dosageMatrix)),]$p_int)

FindGenoDist <- function(training_set,testing_set,index,int_coeff){
  training_geno<- training_set$dosageMatrix[index,]
  geno_dist <- apply(testing_set$dosageMatrix,1,function(x) sum(training_geno == x))
  return(geno_dist)
  
}
# cov_dist <- FindCovDist(training_set,testing_set,1,fit_coeff)
# geno_dist <- FindGenoDist(training_set,testing_set,1,int_coeff)
# 
# mean_rank <- (rank(cov_dist) + rank(geno_dist)) / 2
# 
# training_subset <- which(geno_dist < quantile(geno_dist,probs = 0.05))
# summary(lm(formula = paste0('LDL ~ HMGCR + ',paste(colnames(testing_set$covariates[training_subset,]),collapse = '+')),data = cbind(data.frame(LDL = testing_set$phenotypes[training_subset],HMGCR = testing_set$dosageTarget[training_subset]),testing_set$covariates[training_subset,])))
