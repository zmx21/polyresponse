source('~/MRC_BSU_Internship/Load_Bgen/LoadBgen.R')
source('~/MRC_BSU_Internship/Target_Score/iterative_pruning.R')
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

CalcMarginalEffectGeneScore <- function(path,sample_file_prefix,bgen_file_prefix,phenotype,rsid,eur_only,cov,beta_coeff){
  #Load phenotype information of samples
  print('Loading Phenotypes')
  #Decide what columns to load based on what covariates were specificed
  if(cov!=''){
    cov_names <- unlist(strsplit(x=cov,split = ','))
  }else{
    cov_names <- c()
  }
  #Load phenotype table
  samplePhenoTbl <- data.table::fread(paste0(path,sample_file_prefix,'.csv'),select = c('UKB_genetic_ID','euro',phenotype,cov_names))
  
  #Construct Sample vs Phenotype Table
  phenotypes<- dplyr::select(samplePhenoTbl,phenotype) %>% t() %>% as.vector() 
  if(length(cov_names)>0){
    covariates <- dplyr::select(samplePhenoTbl,cov_names)
    #Remove samples with no phenotype or cov measure.
    samplesToKeep <- !apply(subset(samplePhenoTbl,select=c(phenotype,cov_names)),1,function(x) any(is.na(x)))
    
  }else{
    covariates <- data.frame()
    #Remove samples with no phenotype.
    samplesToKeep <- !apply(subset(samplePhenoTbl,select=c(phenotype)),1,function(x) any(is.na(x)))
    
  }
  #Remove samples with non-european ancestry if specified in argument.
  if(eur_only != 1 & eur_only != 0){
    stop('Please specify eur only')
  }else if(eur_only == 1){
    samplesToKeep <- samplesToKeep & (samplePhenoTbl$euro == 1)
  }
  
  #Remove samples with no phenotyes
  phenotypes <- phenotypes[samplesToKeep]
  if(ncol(covariates) > 0){
    covariates <- as.data.frame(covariates[samplesToKeep,])
  }
  
  #Convert to numeric for non-numeric phenotypes/covariates
  if(!is.numeric(phenotypes)){
    phenotypes <- as.numeric(as.factor(phenotypes))
  }
  if(ncol(covariates) > 0){
    for(cov_name in cov_names){
      curVect <- covariates[,cov_name]
      if(!is.numeric(curVect)){
        covariates[,cov_name] <- as.numeric(as.factor(as.vector(t(curVect))))
      }
    }
  }
  
  dosageVector <- LoadBgen(path,bgen_file_prefix,rsid)
  dosageVector <- dosageVector[,samplesToKeep]
  
  dosageVector <- beta_coeff %*% dosageVector
  fit <- LinearFitGeneScore(as.vector(dosageVector),phenotypes,covariates)
  return(fit)
}
path <-  '~/bsu_scratch/UKB_Data/'
sample_file_prefix <- 'ukbb_metadata'
bgen_file_prefix <- 'ukb_imp_chr#_HRConly'
phenotype <- 'sbp'

distances <- c(5000,seq(10000,200000,10000))
fits <- vector(mode = 'list',length = length(distances))

for(i in 1:length(distances)){
  pruningResult <- IterativePruning('SCNN1D',distances[i],distances[i],5e-6,0.2,20)
  fits[[i]] <- CalcMarginalEffectGeneScore(path,sample_file_prefix,bgen_file_prefix,phenotype,pruningResult$rsid,1,'sex,ages,bmi',pruningResult$coeff)
  write(fits[[i]],ncol=3,file='~/pruning_fits.txt',append=T)
  write(pruningResult$rsid,ncol=length(pruningResult$rsid),file='~/pruning_rsid.txt',append=T)
}
