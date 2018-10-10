library(RcppEigen)
#Function which calculates significance of interaction between two SNPs.
CalcInteractions <- function(dosageSubMatrix,dosageTarget,phenotypes,covariates){
  #Create model matrix, with main and interaction effects of two SNPs. 
  mdlMat <- cbind('Intercept' = rep(1,length(dosageSubMatrix)),'snp' = dosageSubMatrix,'target' = dosageTarget,'int' = dosageSubMatrix * dosageTarget,covariates)
  #Fit linear model and calculate stats
  fit <- RcppEigen::fastLmPure(y = phenotypes,X = mdlMat)
  coeff <- fit$coefficients
  se <- fit$se
  t <- coeff/se
  p <- sapply(t,function(x) 2*pt(abs(x),df=fit$df.residual,lower.tail = F))
  #Form matrix, rows are type of coeff, columns are type of variable.
  fit_result <- rbind(coeff,se,t,p)
  #Return coefficient and significance of interaction term
  return(c(as.vector(fit_result[,-1]),fit$s))
}

RunGxGInteractions <- function(path,sample_file_prefix,bgen_file_prefix,chr,phenotype,targetRS,path_out,eur_only,cov,PC,med,MAF,Info,n_cores){
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
  allRSIds <- FindAllRSIds(chr,as.numeric(MAF),as.numeric(Info)) %>% dplyr::filter(!rsid%in%targetRS)
  allRSIds <- unique(allRSIds$rsid)
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
    
    beta_coeff <- CalcMarginalEffect(path,sample_file_prefix,bgen_file_prefix,phenotype,targetRS,1,cov='sex,ages,bmi',PC=5,med=1,n_cores,F)$coeff
    #Flip alleles such that all are bp lowering
    dosageVectorTemp <- dosageVector
    for(i in 1:nrow(dosageVectorTemp)){
      if(beta_coeff[i] > 0){
        dosageVectorTemp[i,] <- 2-dosageVector[i,]
      }
    }
    dosageVector <- as.vector(abs(beta_coeff) %*% dosageVectorTemp)
  }
  #Run chunks in parallel
  print(paste0('Calculating Interactions: ',length(rsIDChunks),' chunks'))
  allResultsTbl <- pbmclapply(1:length(rsIDChunks),function(i) {
    currentRSIdChunk <- rsIDChunks[[i]]
    dosageMatrix <- LoadBgen(path,bgen_file_prefix,currentRSIdChunk,rep(chr,length(currentRSIdChunk)))
    
    #Keep SNPs other than the target SNP    
    nonTargetSnps <- setdiff(rownames(dosageMatrix),targetRS)
    dosageMatrix <- dosageMatrix[,samplesToKeep]

    #Intialize results vector
    snp_names <- unlist(lapply(c('snp','target','int'),function(x) c(paste0('coeff_',x),paste0('std_err_',x),paste0('t_',x),paste0('p_',x))))
    cov_names <- unlist(lapply(cov_names,function(x) c(paste0('coeff_',x),paste0('std_err_',x),paste0('t_',x),paste0('p_',x))))
    resultsTbl <- read.csv(text = paste(c('rsid',snp_names,cov_names,'RMSE'),collapse = ','))
    resultsTbl[1:length(nonTargetSnps),] <- NA
    a <- Sys.time()
    for(k in 1:length(nonTargetSnps)){
      resultsTbl[k,1] <- nonTargetSnps[k]
      resultsTbl[k,2:ncol(resultsTbl)] <- CalcInteractions(dosageMatrix[nonTargetSnps[k],],dosageTarget,phenotypes,covariates)
    }
    data.table::fwrite(resultsTbl,file = paste0(path_out,'chunk',i,'.txt'),sep = '\t',col.names = T,row.names = F)
  },mc.cores = as.numeric(n_cores),ignore.interactive = T)
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
  bgen_file_prefix <- 'test'
  chr <- '10'
  phenotype = 'sbp'
  targetRS <- 'rs603424'
  target_chr <- '10'
  n_cores <- 1
  eur_only <- 1
  out_suffix <- 'eur_only'
  cov <- 'sex,ages,bmi'
  PC <-  5
  MAF <- 0.05
  info <- 0.5
  med=1
  if(out_suffix == ''){
    path_out <- paste0(path,gsub(',','_',targetRS),'_',phenotype,'/')
    path_out_chr <- paste0(path,gsub(',','_',targetRS),'_',phenotype,'/chr',chr,'/')
  }else{
    path_out <- paste0(path,gsub(',','_',targetRS),'_',phenotype,'_',out_suffix,'/')
    path_out_chr <- paste0(path,gsub(',','_',targetRS),'_',phenotype,'_',out_suffix,'/chr',chr,'/')
  }
  system(paste0('mkdir -p ',path_out))
  system(paste0('chmod a+rwx ',path_out))
  system(paste0('mkdir -p ',path_out_chr))
  system(paste0('chmod a+rwx ',path_out_chr))

}else if (length(args)!= 14){
  stop("You need to supply:\n",
       "# 1: Input Path\n",
       '# 2: Sample File Prefix\n',
       "# 3: Bgen File Prefix\n",
       "# 4: Chromosome to Analyze\n",
       "# 5: Phenotype Identifier (column name in phenotype file) \n",
       "# 6: rsid of Target SNP\n",
       "# 7: Suffix of output directory\n",
       "# 8: EUR only\n",
       "# 9: Covariates (seperated by comma)\n",
       "# 10: Number of PC to include\n",
       '# 11: Medication: (add 10 to bp for those taking med)',
       "# 12: MAF filter\n",
       "# 13: Info filter\n",
       '# 14: Number of Cores\n',
       "Exiting...", call.=FALSE)
  
}else{
  print('All Arguments Supplied')
  str <- c("# 1: Input Path:",
           '# 2: Sample File Prefix:',
           "# 3: Bgen File Prefix:",
           "# 4: Chromosome to Analyze:",
           "# 5: Phenotype Identifier (column name in phenotype file):",
           "# 6: rsid of Target SNP:",
           "# 7: Suffix of output directory:",
           "# 8: EUR only:",
           "# 9: Covariates (seperated by comma):",
           "# 10: Number of PC to include",
           '# 11: Medication: (add 10 to bp for those taking med)',
           "# 12: MAF cutoff",
           "# 13: Info cutoff",
           '# 14: Number of Cores:')
  variables <- c('path','sample_file_prefix','bgen_file_prefix','chr','phenotype','targetRS','out_suffix','eur_only','cov','PC','med','MAF','info','n_cores')
  for(i in 1:length(args)){
    eval(parse(text=paste0(variables[i],'=',"'",args[[i]],"'")))
    print(paste0(str[i],"   ",args[i]))
  }
  if(out_suffix == ''){
    path_out <- paste0(path,gsub(',','_',targetRS),'_',phenotype,'/')
    path_out_chr <- paste0(path,gsub(',','_',targetRS),'_',phenotype,'/chr',chr,'/')
  }else{
    path_out <- paste0(path,gsub(',','_',targetRS),'_',phenotype,'_',out_suffix,'/')
    path_out_chr <- paste0(path,gsub(',','_',targetRS),'_',phenotype,'_',out_suffix,'/chr',chr,'/')
  }
  system(paste0('mkdir -p ',path_out))
  system(paste0('chmod a+rwx ',path_out))
  system(paste0('mkdir -p ',path_out_chr))
  system(paste0('chmod a+rwx ',path_out_chr))
}
RunGxGInteractions(path,sample_file_prefix,bgen_file_prefix,chr,phenotype,targetRS,path_out_chr,eur_only,
                   cov,PC,med,MAF,info,n_cores)
