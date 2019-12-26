source('~/MRC_BSU_Internship_LDL/Load_Bgen/LoadBgen.R')
source('~/MRC_BSU_Internship_LDL/Load_Phenotype/Load_Phenotype.R')

library(RcppEigen)
library(pbmcapply)
library(dplyr)
LinearFit <- function(dosageSubMatrix,phenotypes,covariates){  
#Create model matrix, with main and interaction effects of two SNPs. 
  if(ncol(covariates) > 0){
    mdlMat <- cbind('Intercept' = rep(1,length(dosageSubMatrix)),'snp' = dosageSubMatrix,covariates)
  }else{
    mdlMat <- cbind('Intercept' = rep(1,length(dosageSubMatrix)),'snp' = dosageSubMatrix)
  }
  
  fit <- RcppEigen::fastLm(y = phenotypes,X = mdlMat)
  #Return coefficient and significance of interaction term
  results <- summary(fit)$coefficients[2,c(1,2,4)]
  names(results) <- c('coeff','se','p')
  return(list(results = results,fit=summary(fit)))
}

CalcSnpPhenoAssociation <- function(path,sample_file_prefix,bgen_file_prefix,phenotype,rsid,eur_only,cov,PC,med,n_cores,verbose=F){
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
  
  print('Calculating Fit')
  if(verbose == T){
    results <- vector(mode = 'list',length = length(rsid))
    for(i in 1:length(rsid)){
      dosageVector <- LoadBgen(path,bgen_file_prefix,rsid[i])
      dosageVector <- dosageVector[,samplesToKeep]
      results[[i]] <- LinearFit(dosageVector,phenotypes,covariates)$fit
      names(results) <- rsid
      return(results)
    }
  }else{
    #Calc marginal effects
    results <- lapply(1:length(rsid),function(i) {
      fit = tryCatch({
        dosageVector <- LoadBgen(path,bgen_file_prefix,rsid[i])
        dosageVector <- dosageVector[,samplesToKeep]
        LinearFit(dosageVector,phenotypes,covariates)$results
      }, warning = function(w) {
        return(c(coeff=NA,p=NA))
      }, error = function(e) {
        return(c(coeff=NA,p=NA))
      })
      return(fit)
    })%>% {do.call(rbind,.)}
    results <- as.data.frame(results)
    results <- cbind(data.frame(rsid=rsid,stringsAsFactors = F),results)
    return(results)
  }
}
