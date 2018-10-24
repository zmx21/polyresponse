source('~/MRC_BSU_Internship/Load_Bgen/LoadBgen.R')
# source('~/MRC_BSU_Internship/Target_Score/iterative_pruning.R')
source('~/MRC_BSU_Internship/Load_Phenotype/Load_Phenotype.R')

library(RcppEigen)
library(pbmcapply)
library(dplyr)
LinearFitGeneScore <- function(dosageSubMatrix,phenotypes,covariates){
  #Create model matrix, with main and interaction effects of two SNPs. 
  if(ncol(covariates) > 0){
    mdlMat <- cbind('Intercept' = rep(1,length(dosageSubMatrix)),'gene_score' = dosageSubMatrix,covariates)
  }else{
    mdlMat <- cbind('Intercept' = rep(1,length(dosageSubMatrix)),'gene_score' = dosageSubMatrix)
  }
  
  fit <- RcppEigen::fastLm(y = phenotypes,X = mdlMat)
  #Return coefficient and significance of interaction term
  results <- c(summary(fit)$coefficients[2,c(1,4)],summary(fit)$r.squared)
  names(results) <- c('coeff','p','rsq')
  return(results)
}

CalcMarginalEffectGeneScore <- function(path,sample_file_prefix,bgen_file_prefix,phenotype,rsid,eur_only,cov,beta_coeff,med=T,PC=5){
  #Decide what columns to load based on what covariates were specificed
  if(cov!=''){
    cov_names <- unlist(strsplit(x=cov,split = ','))
  }else{
    cov_names <- c()
  }
  #Add PC to covariates, if specified
  if(PC > 0){
    cov_names <- c(cov_names,sapply(1:PC,function(x) paste0('PC',x)))
  }
  
  #Load Phenotype and Covariates
  print('Loading Phenotypes')
  phenotypesAndCov <- LoadPhenotype(path,sample_file_prefix,phenotype,cov_names,eur_only,med)
  phenotypes <- phenotypesAndCov$phenotypes
  covariates <- phenotypesAndCov$covariates
  samplesToKeep <- phenotypesAndCov$samplesToKeep
  
  dosageVector <- LoadBgen(path,bgen_file_prefix,rsid)
  dosageVector <- dosageVector[,samplesToKeep]
  
  dosageVector <- beta_coeff %*% dosageVector
  fit <- LinearFitGeneScore(as.vector(dosageVector),phenotypes,covariates)
  return(fit)
}
# path <-  '~/bsu_scratch/UKB_Data/'
# sample_file_prefix <- 'ukbb_metadata'
# bgen_file_prefix <- 'ukb_imp_chr#_HRConly'
# phenotype <- 'sbp'

target_names <- c('SCNN1A','SCNN1B','SCNN1G','SCNN1D','ACE','CACNA1D','CACNA1S','CACNA1H','CACNA2D1','CACNA2D2','AGTR1','MME','ADRA2B','GPAT2','PDE5A')
target_sbp <- lapply(target_names,function(x){
  result <- IterativePruning(x,'sbp',10000,10000,5e-6,0.3,12)
  cbind(data.frame(gene_name=rep(x,nrow(result$coeff)),stringsAsFactors = F),result$coeff)})
data.table::fwrite(do.call(rbind,target_sbp),row.names = F,file = '~/bsu_scratch/target_sbp_ukbld.txt')
target_dbp <- lapply(target_names,function(x){
  result <- IterativePruning(x,'dbp',10000,10000,5e-6,0.3,12)
  cbind(data.frame(gene_name=rep(x,nrow(result$coeff)),stringsAsFactors = F),result$coeff)})
data.table::fwrite(do.call(rbind,target_dbp),row.names = F,file = '~/bsu_scratch/target_dbp_ukbld.txt')


# SCNN1D_sbp <- IterativePruning('SCNN1D','sbp',10000,10000,5e-6,0.3,20)
# ACE_sbp <- IterativePruning('ACE','sbp',10000,10000,5e-6,0.3,20)
# CACNA2D2_sbp<- IterativePruning('CACNA2D2','sbp',10000,10000,5e-6,0.3,20)
# MME_sbp <- IterativePruning('MME','sbp',10000,10000,5e-6,0.3,20)
# ADRA2B_sbp <- IterativePruning('ADRA2B','sbp',10000,10000,5e-6,0.3,20)
# GPAT2_sbp <- IterativePruning('GPAT2','sbp',10000,10000,5e-6,0.3,20)
# PDE5A_sbp <- IterativePruning('PDE5A','sbp',10000,10000,5e-6,0.3,20)
# 
# SCNN1D_dbp <- IterativePruning('SCNN1D','dbp',10000,10000,5e-6,0.3,20)
# ACE_dbp <- IterativePruning('ACE','dbp',10000,10000,5e-6,0.3,20)
# CACNA2D2_dbp<- IterativePruning('CACNA2D2','dbp',10000,10000,5e-6,0.3,20)
# MME_dbp <- IterativePruning('MME','dbp',10000,10000,5e-6,0.3,20)
# ADRA2B_dbp <- IterativePruning('ADRA2B','dbp',10000,10000,5e-6,0.3,20)
# GPAT2_dbp <- IterativePruning('GPAT2','dbp',10000,10000,5e-6,0.3,20)
# PDE5A_dbp <- IterativePruning('PDE5A','dbp',10000,10000,5e-6,0.3,20)

# distances <- c(5000,seq(10000,200000,10000))
# fits <- vector(mode = 'list',length = length(distances))
# 
# for(i in 1:length(distances)){
#   pruningResult <- IterativePruning('SCNN1D',distances[i],distances[i],5e-6,0.2,20)
#   fits[[i]] <- CalcMarginalEffectGeneScore(path,sample_file_prefix,bgen_file_prefix,phenotype,pruningResult$rsid,1,'sex,ages,bmi',pruningResult$coeff)
#   write(fits[[i]],ncol=3,file='~/pruning_fits.txt',append=T)
#   write(pruningResult$rsid,ncol=length(pruningResult$rsid),file='~/pruning_rsid.txt',append=T)
# }
