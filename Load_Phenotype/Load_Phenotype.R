#Load phenotype table, according to filters and covariates specified.
library(dplyr)
LoadPhenotype <- function(path,sample_file_prefix,phenotype,cov_names,eur_only,med,verbose=F,excl=c()){
  #Load phenotype table
  samplePhenoTbl <- data.table::fread(paste0(path,sample_file_prefix,'.csv'),
                                      select = c('UKB_genetic_ID','euro','lipdbin',phenotype,cov_names))
  if(verbose){
    print(paste0('Raw Samples ',nrow(samplePhenoTbl)))
    print(paste0('Euro Samples',nrow(dplyr::filter(samplePhenoTbl,euro==1))))
  }
  
  if(med != 1 & med != 0){
    stop('Please specify med')
  }
  #Multiply 1.25 to LDL for medicine takers (medication reduces LDL by 20%).
  if(med == 1){
    medicineTakers <- samplePhenoTbl$lipdbin == 'Current'
    samplePhenoTbl$LDLdirect[medicineTakers] <- samplePhenoTbl$LDLdirect[medicineTakers] * 1.25
  }
  
  #Construct Sample vs Phenotype Table
  phenotypes<- dplyr::select(samplePhenoTbl,phenotype) %>% t() %>% as.vector() 
  if(length(cov_names)>0){
    covariates <- dplyr::select(samplePhenoTbl,cov_names)
    #Remove samples with no phenotype or cov measure.
    #Keep missing rows for cov not included in analysis
    cov_names_analysis <- setdiff(cov_names,excl)
    samplesToKeep <- !apply(subset(samplePhenoTbl,select=c(phenotype,cov_names_analysis)),1,function(x) any(is.na(x)))
    
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
  
  #Convert to binary for categorical variables with 2 category
  if(!is.numeric(phenotypes)){
    phenotypes <- as.numeric(as.factor(phenotypes)) - 1
  }
  if(ncol(covariates) > 0){
    for(cov_name in cov_names){
      curVect <- covariates[,cov_name]
      if(!is.numeric(curVect)){
        covariates[,cov_name] <- as.numeric(as.factor(as.vector(t(curVect)))) - 1
      }
    }
  }
  if(verbose){
    print(paste0('Non Missing Covariates: ',length(phenotypes)))
  }
  
  return(list(phenotypes=phenotypes,covariates=covariates,samplesToKeep=samplesToKeep,lipdbin = samplePhenoTbl$lipdbin[samplesToKeep]))
}

