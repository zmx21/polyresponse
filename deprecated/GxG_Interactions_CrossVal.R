library(RcppEigen)
#Parse Fit result
ParseFit <- function(fit){
  coeff <- fit$coefficients
  se <- fit$se
  t <- coeff/se
  p <- sapply(t,function(x) 2*pt(abs(x),df=fit$df.residual,lower.tail = F))
  #Form matrix, rows are type of coeff, columns are type of variable.
  fit_result <- rbind(coeff,se,t,p)
  return(fit_result)
}

#Function which calculates significance of interaction between two SNPs.
MultipleRegression <- function(dosageSubMatrix,dosageTarget,phenotypes,covariates){
  #Create model matrix, with main and interaction effects of two SNPs. 
  mdlMatInt <- cbind('Intercept' = rep(1,length(dosageSubMatrix)),'snp' = dosageSubMatrix,'target' = dosageTarget,'int' = dosageSubMatrix * dosageTarget,covariates)
  #Create model matrix, with main effects of two SNPs. 
  mdlMatNoInt <- cbind('Intercept' = rep(1,length(dosageSubMatrix)),'snp' = dosageSubMatrix,'target' = dosageTarget,covariates)
  
  #Fit linear model and calculate stats
  fitInt <- RcppEigen::fastLmPure(y = phenotypes,X = mdlMatInt)
  fitNoInt <- RcppEigen::fastLmPure(y = phenotypes,X = mdlMatNoInt)
  
  intResults <- ParseFit(fitInt)
  noIntResults <- ParseFit(fitNoInt)
}

RunGxGInteractionsCrossVal <- function(path,sample_file_prefix,bgen_file_prefix,phenotype,targetRS,interaction_results_path,interaction_thresh,path_out,eur_only,cov,PC,med,k_fold,n_cores){
  chunkSize=50
  
  library(dplyr)
  library(parallel)
  library(data.table)
  library(pbmcapply)
  source('~/MRC_BSU_Internship/Load_Bgen/LoadBgen.R')
  source('~/MRC_BSU_Internship/Load_Phenotype/Load_Phenotype.R')
  source('~/MRC_BSU_Internship/SNP_Marginal_Effect/CalcMarginalEffect.R')
  #Find all rsids, and generate chunks to read.
  print('Loading rsID')
  interactionResults <- system(paste0("awk '{if($13 < ",interaction_thresh,"|| $13==\"p_int\") {print$1\"\\t\"$13} }' ",interaction_results_path),intern = T)
  interactionResults <- data.table::fread(text = interactionResults,header = T,sep = '\t')
  interactionResults$p_int <- as.numeric(interactionResults$p_int)
  allRSIds <- interactionResults %>% dplyr::filter(p_int != 0) %>% {unique(.$rsid)}
  rsIDChunks <- split(allRSIds,seq(length(allRSIds)-1)%/%chunkSize)
  #Load phenotype information of samples
  print('Loading Phenotypes')
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
  phenotypesAndCov <- LoadPhenotype(path,sample_file_prefix,phenotype,cov_names,eur_only,med)
  phenotypes <- phenotypesAndCov$phenotypes
  covariates <- as.matrix(phenotypesAndCov$covariates)
  samplesToKeep <- phenotypesAndCov$samplesToKeep
  
  
  #Split rsid based on comma
  targetRS <- unlist(strsplit(targetRS,split = ','))
  
  #Genotype info of target gene
  if(length(targetRS) < 2){
    dosageTarget <- LoadBgen(path,bgen_file_prefix,targetRS)
    dosageTarget <- dosageTarget[,samplesToKeep]
  }else{
    #Get marginal effects of rsid
    dosageVector <- LoadBgen(path,bgen_file_prefix,targetRS)
    dosageVector <- dosageVector[,samplesToKeep]
    
    beta_coeff <- CalcMarginalEffect(path,sample_file_prefix,bgen_file_prefix,phenotype,targetRS,as.numeric(eur_only),cov=cov,PC=as.numeric(PC),med=as.numeric(med),as.numeric(n_cores),F)$coeff
    #Flip alleles such that all are bp lowering
    dosageTarget <- dosageVector
    for(i in 1:nrow(dosageTarget)){
      if(beta_coeff[i] > 0){
        dosageTarget[i,] <- 2-dosageVector[i,]
      }
    }
    dosageTarget <- as.vector(abs(beta_coeff) %*% dosageTarget)
  }
  
  #Run chunks in parallel
  # print(paste0('Calculating Interactions: ',length(rsIDChunks),' chunks'))
  # allResultsTbl <- pbmclapply(1:length(rsIDChunks),function(i) {
    i <- 1
    currentRSIdChunk <- rsIDChunks[[i]]
    dosageMatrix <- LoadBgen(path,bgen_file_prefix,currentRSIdChunk)
    
    #Keep SNPs other than the target SNP    
    nonTargetSnps <- setdiff(rownames(dosageMatrix),targetRS)
    dosageMatrix <- dosageMatrix[,samplesToKeep]

    #Intialize results vector
    resultsTbl <- data.frame(rsid=character(),mean_coeff_int=numeric(),sd_coeff_int=numeric(),
                             mean_p_int=numeric(),sd_p_int=numeric(),
                             mean_RMSE_int=numeric(),sd_RMSE_int=numeric(),
                             mean_RMSE_no_int=numeric(),sd_RMSE_no_int=numeric())
    resultsTbl[1:length(nonTargetSnps),] <- NA
    for(k in 1:length(nonTargetSnps)){
      resultsTbl[k,1] <- nonTargetSnps[k]
      resultsTbl[k,2:ncol(resultsTbl)] <- MultipleRegression(dosageMatrix[nonTargetSnps[k],],dosageTarget,phenotypes,covariates)
    }
  #   data.table::fwrite(resultsTbl,file = paste0(path_out,'chunk',i,'.txt'),sep = '\t',col.names = T,row.names = F)
  # },mc.cores = as.numeric(n_cores),ignore.interactive = T)
}
##First read in the arguments listed at the command line
args=(commandArgs(TRUE))
##args is now a list of character vectors
## First check to see if arguments are passed.
## Then cycle through each element of the list and evaluate the expressions.
if(length(args)==0){
  print("No arguments supplied.")
  #supply default values
  path <-  '/mrc-bsu/scratch/zmx21/UKB_Data/'
  sample_file_prefix <- 'ukbb_metadata_with_PC'
  bgen_file_prefix <- 'ukb_imp_chr#_HRConly'
  phenotype = 'sbp'
  targetRS <- 'rs1262894'
  interaction_results_path <- paste0('~/bsu_scratch/UKB_Data/',targetRS,'_',phenotype,'_fullcorrection/all_chr.txt')
  interaction_thresh <- 0.01
  eur_only <- 1
  out_suffix <- 'test'
  cov <- 'sex,ages,bmi'
  PC <-  5
  MAF <- 0.05
  info <- 0.5
  med <- 1
  k_fold <- 5
  n_cores <- 1
  if(out_suffix == ''){
    path_out <- paste0(path,gsub(',','_',targetRS),'_',phenotype,'_CV','/')
  }else{
    path_out <- paste0(path,gsub(',','_',targetRS),'_',phenotype,'_CV_',out_suffix,'/')
  }
  system(paste0('mkdir -p ',path_out))
  system(paste0('chmod a+rwx ',path_out))
}else if (length(args)!= 13){
  stop("You need to supply:\n",
       "# 1: Input Path:\n",
       '# 2: Sample File Prefix:\n',
       "# 3: Bgen File Prefix:\n",
       "# 4: Phenotype Identifier (column name in phenotype file):\n",
       "# 5: rsid of Target SNP:\n",
       "# 6: Interaction threshold to run CV:\n",
       "# 7: Suffix of output directory:\n",
       "# 8: EUR only:\n",
       "# 9: Covariates (seperated by comma):\n",
       "# 10: Number of PC to include\n",
       '# 11: Medication: (add 10 to bp for those taking med):\n',
       '# 12: Number of folds for CV:\n',
       '# 13: Number of Cores:\n',
       "Exiting...", call.=FALSE)
  
}else{
  print('All Arguments Supplied')
  str <- c("# 1: Input Path:",
           '# 2: Sample File Prefix:',
           "# 3: Bgen File Prefix:",
           "# 4: Phenotype Identifier (column name in phenotype file):",
           "# 5: rsid of Target SNP:",
           "# 6: Interaction threshold to run CV:",
           "# 7: Suffix of output directory:",
           "# 8: EUR only:",
           "# 9: Covariates (seperated by comma):",
           "# 10: Number of PC to include:",
           '# 11: Medication: (add 10 to bp for those taking med):',
           '# 12: Number of folds for CV:',
           '# 13: Number of Cores:')
  variables <- c('path','sample_file_prefix','bgen_file_prefix','phenotype','targetRS','interaction_thresh','eur_only','cov','PC','med','k_fold','n_cores')
  for(i in 1:length(args)){
    eval(parse(text=paste0(variables[i],'=',"'",args[[i]],"'")))
    print(paste0(str[i],"   ",args[i]))
  }
  interaction_results_path <- paste0('~/bsu_scratch/UKB_Data/',targetRS,'_',phenotype,'_fullcorrection/all_chr.txt')
  if(out_suffix == ''){
    path_out <- paste0(path,gsub(',','_',targetRS),'_',phenotype,'_CV','/')
  }else{
    path_out <- paste0(path,gsub(',','_',targetRS),'_',phenotype,'_CV_',out_suffix,'/')
  }
  system(paste0('mkdir -p ',path_out))
  system(paste0('chmod a+rwx ',path_out))
}
RunGxGInteractionsCrossVal(path,sample_file_prefix,bgen_file_prefix,phenotype,targetRS,interaction_results_path,interaction_thresh,path_out,eur_only,cov,PC,med,k_fold,n_cores)
