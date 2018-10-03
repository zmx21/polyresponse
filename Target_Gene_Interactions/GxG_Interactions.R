library(RcppEigen)
#Function which calculates significance of interaction between two SNPs.
CalcInteractions <- function(dosageSubMatrix,dosageTarget,phenotypes,covariates){
  #Create model matrix, with main and interaction effects of two SNPs. 
  mdlMat <- cbind('Intercept' = rep(1,length(dosageSubMatrix)),'snp' = dosageSubMatrix,'target' = dosageTarget,'int' = dosageSubMatrix * dosageTarget,covariates)
  fit <- RcppEigen::fastLm(y = phenotypes,X = mdlMat)
  #Return coefficient and significance of interaction term
  return(c(as.vector(t(summary(fit)$coefficients[-1,])),summary(fit)$r.squared))
}

RunGxGInteractions <- function(path,sample_file_prefix,bgen_file_prefix,chr,phenotype,targetRS,path_out,eur_only,cov,n_cores){
  chunkSize=50
  
  library(dplyr)
  library(parallel)
  library(data.table)
  library(pbmcapply)
  source('/home/zmx21/MRC_BSU_Internship/Load_Bgen/LoadBgen.R')
  
  #Genotype info of target gene
  dosageTarget <- LoadBgen(path,bgen_file_prefix,targetRS)

  #Find all rsids, and generate chunks to read.
  print('Loading rsID')
  allRSIds <- FindAllRSIds(chr) %>% dplyr::filter(rsid!=targetRS)
  allRSIds <- unique(allRSIds$rsid)
  rsIDChunks <- split(allRSIds,seq(length(allRSIds)-1)%/%chunkSize)
  #Load phenotype information of samples
  print('Loading Phenotypes')
  #Decide what columns to load based on what covariates were specificed
  cov_names <- unlist(strsplit(x=cov,split = ','))
  #Load phenotype table
  samplePhenoTbl <- data.table::fread(paste0(path,sample_file_prefix,'.csv'),select = c('UKB_genetic_ID','euro',phenotype,cov_names))
  
  #Construct Sample vs Phenotype Table
  phenotypes<- dplyr::select(samplePhenoTbl,phenotype) %>% t() %>% as.vector() 
  covariates <- dplyr::select(samplePhenoTbl,cov_names)
  #Remove samples with no phenotype or cov measure, and non-european ancestry if specified in argument.
  samplesToKeep <- !apply(subset(samplePhenoTbl,select=c(phenotype,cov_names)),1,function(x) any(is.na(x)))
  if(eur_only != 1 & eur_only != 0){
    stop('Please specify eur only')
  }else if(eur_only == 1){
    samplesToKeep <- samplesToKeep & (samplePhenoTbl$euro == 1)
  }
  dosageTarget <- dosageTarget[,samplesToKeep]
  phenotypes <- phenotypes[samplesToKeep]
  covariates <- as.data.frame(covariates[samplesToKeep,])
  
  #Convert to numeric for non-numeric phenotypes/covariates
  if(!is.numeric(phenotypes)){
    phenotypes <- as.numeric(as.factor(phenotypes))
  }
  for(cov_name in cov_names){
    curVect <- covariates[,cov_name]
    if(!is.numeric(curVect)){
      covariates[,cov_name] <- as.numeric(as.factor(as.vector(t(curVect))))
    }
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
    # resultsTbl <- data.frame('rsid'=rep('',length(nonTargetSnps)),'p' = rep(NA,length(nonTargetSnps)), 'coeff'=rep(NA,length(nonTargetSnps)),'rsq'=rep(NA,length(nonTargetSnps)),stringsAsFactors = F)
    snp_names <- unlist(lapply(c('snp','target','int'),function(x) c(paste0('coeff_',x),paste0('std_err_',x),paste0('t_',x),paste0('p_',x))))
    cov_names <- unlist(lapply(cov_names,function(x) c(paste0('coeff_',x),paste0('std_err_',x),paste0('t_',x),paste0('p_',x))))
    resultsTbl <- read.csv(text = paste(c('rsid',snp_names,cov_names,'rsq'),collapse = ','))
    resultsTbl[1:length(nonTargetSnps),] <- NA
    for(k in 1:length(nonTargetSnps)){
      resultsTbl[k,1] <- nonTargetSnps[k]
      resultsTbl[k,2:ncol(resultsTbl)] <- CalcInteractions(dosageMatrix[nonTargetSnps[k],],dosageTarget,phenotypes,covariates)
    }
    data.table::fwrite(resultsTbl,file = paste0(path_out,'chunk',i,'.txt'),sep = '\t',col.names = F,row.names = F)
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
  sample_file_prefix <- 'ukbb_metadata'
  bgen_file_prefix <- 'ukb_imp_chr#_HRConly'
  chr <- '2'
  phenotype = 'sbp'
  targetRS <- 'rs6602047'
  target_chr <- '10'
  n_cores=1
  eur_only=1
  out_suffix <- 'eur_only'
  cov='sex,ages,bmi'
  if(out_suffix == ''){
    path_out <- paste0(path,targetRS,'_',phenotype,'/')
    path_out_chr <- paste0(path,targetRS,'_',phenotype,'/chr',chr,'/')
  }else{
    path_out <- paste0(path,targetRS,'_',phenotype,'_',out_suffix,'/')
    path_out_chr <- paste0(path,targetRS,'_',phenotype,'_',out_suffix,'/chr',chr,'/')
  }
  system(paste0('mkdir -p ',path_out))
  system(paste0('chmod a+rwx ',path_out))
  system(paste0('mkdir -p ',path_out_chr))
  system(paste0('chmod a+rwx ',path_out_chr))

}else if (length(args)!=10){
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
       '# 10: Number of Cores\n',
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
           '# 10: Number of Cores:')
  variables <- c('path','sample_file_prefix','bgen_file_prefix','chr','phenotype','targetRS','out_suffix','eur_only','cov','n_cores')
  for(i in 1:length(args)){
    eval(parse(text=paste0(variables[i],'=',"'",args[[i]],"'")))
    print(paste0(str[i],"   ",args[i]))
  }
  if(out_suffix == ''){
    path_out <- paste0(path,targetRS,'_',phenotype,'/')
    path_out_chr <- paste0(path,targetRS,'_',phenotype,'/chr',chr,'/')
  }else{
    path_out <- paste0(path,targetRS,'_',phenotype,'_',out_suffix,'/')
    path_out_chr <- paste0(path,targetRS,'_',phenotype,'_',out_suffix,'/chr',chr,'/')
  }
  system(paste0('mkdir -p ',path_out))
  system(paste0('chmod a+rwx ',path_out))
  system(paste0('mkdir -p ',path_out_chr))
  system(paste0('chmod a+rwx ',path_out_chr))
}
RunGxGInteractions(path,sample_file_prefix,bgen_file_prefix,chr,phenotype,targetRS,path_out_chr,eur_only,cov,n_cores)
