#Load phenotype table, according to filters and covariates specified.
LoadPhenotype <- function(path,sample_file_prefix,phenotype,cov_names,eur_only,med){
  #Load phenotype table
  samplePhenoTbl <- data.table::fread(paste0(path,sample_file_prefix,'.csv'),
                                      select = c('UKB_genetic_ID','euro','hypdbin',phenotype,cov_names))
  #Add 10 to blood pressure of medicine takers.
  if(eur_only != 1 & eur_only != 0){
    stop('Please specify med')
  }else if(med == 1 & phenotype == 'sbp'){
    medicineTakers <- samplePhenoTbl$hypdbin == 'Current'
    samplePhenoTbl$sbp[medicineTakers] <- samplePhenoTbl$sbp[medicineTakers] + 10
  }else if(med == 1 & phenotype == 'dbp'){
    medicineTakers <- samplePhenoTbl$hypdbin == 'Current'
    samplePhenoTbl$dbp[medicineTakers] <- samplePhenoTbl$dbp[medicineTakers] + 10
  }
  
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
  
  return(list(phenotypes=phenotypes,covariates=covariates,samplesToKeep=samplesToKeep))
}