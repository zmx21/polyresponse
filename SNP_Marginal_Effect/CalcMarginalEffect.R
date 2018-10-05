source('/home/zmx21/MRC_BSU_Internship/Load_Bgen/LoadBgen.R')
source('~/MRC_BSU_Internship/Load_Phenotype/Load_Phenotype.R')

library(RcppEigen)
library(pbmcapply)
library(dplyr)
LinearFit <- function(dosageSubMatrix,phenotypes,covariates,verbose=T){
  #Create model matrix, with main and interaction effects of two SNPs. 
  if(ncol(covariates) > 0){
    mdlMat <- cbind('Intercept' = rep(1,length(dosageSubMatrix)),'snp' = dosageSubMatrix,covariates)
  }else{
    mdlMat <- cbind('Intercept' = rep(1,length(dosageSubMatrix)),'snp' = dosageSubMatrix)
  }
  
  fit <- RcppEigen::fastLm(y = phenotypes,X = mdlMat)
  if(verbose){
    summary(fit)
  }
  #Return coefficient and significance of interaction term
  results <- summary(fit)$coefficients[2,c(1,4)]
  names(results) <- c('coeff','p')
  return(results)
}

CalcMarginalEffect <- function(path,sample_file_prefix,bgen_file_prefix,phenotype,rsid,eur_only,cov,PC,med,n_cores,verbose){
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
  
  #Calc marginal effects
  results <- pbmclapply(1:length(rsid),function(i) {
    fit = tryCatch({
      dosageVector <- LoadBgen(path,bgen_file_prefix,rsid[i])
      dosageVector <- dosageVector[,samplesToKeep]
      LinearFit(dosageVector,phenotypes,covariates)
    }, warning = function(w) {
      return(c(coeff=NA,p=NA))
    }, error = function(e) {
      return(c(coeff=NA,p=NA))
    })
    return(fit)
  },mc.cores=n_cores) %>% {do.call(rbind,.)}
  results <- as.data.frame(results)
  results <- cbind(data.frame(rsid=rsid,stringsAsFactors = F),results)
  return(results)
}

# path <-  '/mrc-bsu/scratch/zmx21/UKB_Data/'
# sample_file_prefix <- 'ukbb_metadata'
# bgen_file_prefix <- 'ukb_imp_chr#_HRConly'
# phenotype <- 'sbp'
# 
# with_cov <- CalcMarginalEffect(path,sample_file_prefix,bgen_file_prefix,phenotype,'rs1262894',1,'sex,ages,bmi',1)
# no_cov <- CalcMarginalEffect(path,sample_file_prefix,bgen_file_prefix,phenotype,'rs1262894',1,'',1)
