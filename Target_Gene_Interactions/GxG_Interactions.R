####################################################################################
#Main function to run genomewide GxG interactions.
#Input: As listed in  message below 
#Output: Text file with GxG interaction information, written to a text file in specified path
####################################################################################

library(RcppEigen)
library(dplyr)
library(parallel)
library(data.table)
library(pbmcapply)
library(dbplyr)

source('../Target_Gene_Interactions/CalcInteractions.R')
source('../Load_Bgen/LoadBgen.R')
source('../Load_Phenotype/Load_Phenotype.R')
source('../SNP_Phenotype_Association/CalcSnpPhenoAssociation.R')

#Main function, which calculates GxG interactions with target SNP, in chunks.
RunGxGInteractions <- function(path,sample_file_prefix,bgen_file_prefix,chr,phenotype,targetRS,path_out,eur_only,cov,PC,med,MAF,Info,n_cores,chunks,training_set){
  chunkSize=50
  
  #Find all rsids, and generate chunks to read.
  print('Loading rsID')
  allRSIds <- data.table::fread(paste0(path,'ukb_imp_chr',chr,'_HRConly_samplefilt_snpstats.txt'),skip=10,header = T) %>% dplyr::filter(minor_allele_frequency > as.numeric(MAF)) %>% dplyr::filter(info > as.numeric(Info)) %>% dplyr::filter(!rsid%in%targetRS) #%>% dplyr::filter(rsid == 'rs162724') 
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

  #Load training set index,subset if specified
  if(training_set != 'NULL'){
     trainingSet <- readRDS(file = training_set)
     phenotypes <- phenotypes[trainingSet]
     covariates <- covariates[trainingSet,]
  }else{
     trainingSet <- 1:length(phenotypes)
  }
  
  #Calculate gene score if more than 1 target SNP listed, otherwise load single SNP dosage vector.
  if(length(targetRS) < 2){
    #Load single SNP dosage vector
    dosageTarget <- LoadBgen(path,bgen_file_prefix,targetRS)
    #keep samples without missing information
    dosageTarget <- dosageTarget[,samplesToKeep]
    #keep only samples in training set
    dosageTarget <- dosageTarget[trainingSet]
    
  }else{
    #Get dosage vectors of all SNPs
    dosageVector <- LoadBgen(path,bgen_file_prefix,targetRS)
    #keep samples without missing information
    dosageVector <- dosageVector[,samplesToKeep]
    #keep only samples in training set
    dosageVector <- dosageVector[,trainingSet]
    
    #Calculate main effects of each individual SNP
    beta_coeff <- CalcSnpPhenoAssociation(path,sample_file_prefix,bgen_file_prefix,phenotype,targetRS,
                                          as.numeric(eur_only),cov=cov,PC=as.numeric(PC),
                                          med=as.numeric(med),as.numeric(n_cores),F)$coeff
    #Flip alleles such that all are LDL lowering
    dosageTarget <- dosageVector
    for(i in 1:nrow(dosageTarget)){
      if(beta_coeff[i] > 0){
        dosageTarget[i,] <- 2-dosageVector[i,]
      }
    }
    dosageTarget <- as.vector(abs(beta_coeff) %*% dosageTarget)
  }
  
  #Parse chunks argument
  chunks <- unlist(strsplit(chunks,','))
  if(chunks[2]=='end'){
   chunks <- as.numeric(chunks[1]):length(rsIDChunks) 
  }else{
    chunks <- as.numeric(chunks)
    chunks <- chunks[1]:chunks[2]
  }
  
  #Run chunks in parallel
  print(paste0('Calculating Interactions: ',length(chunks),' chunks'))
  allResultsTbl <- pbmclapply(chunks,function(i){ #lapply(chunks, function(i){  
    # i <- 1
    currentRSIdChunk <- rsIDChunks[[i]]
    dosageMatrix <- LoadBgen(path,bgen_file_prefix,currentRSIdChunk,rep(chr,length(currentRSIdChunk)))
    
    #Keep SNPs other than the target SNP    
    nonTargetSnps <- setdiff(rownames(dosageMatrix),targetRS)
    #remove duplicated SNPs
    loadedSnps <- rownames(dosageMatrix)
    dupSnps <- loadedSnps[which(duplicated(loadedSnps))]
    snpsToKeep <- setdiff(nonTargetSnps,dupSnps)
    #Remove target snps and dup snps
    dosageMatrix <- dosageMatrix[loadedSnps %in% snpsToKeep,]
    #Keep samples without missing info
    dosageMatrix <- dosageMatrix[,samplesToKeep]
    #Keep samples in training set
    dosageMatrix <- dosageMatrix[,trainingSet]
    
    #Intialize results vector
    snp_names <- unlist(lapply(c('snp','target','int'),function(x) c(paste0('coeff_',x),paste0('std_err_',x),paste0('t_',x),paste0('p_',x))))
    cov_names <- unlist(lapply(cov_names,function(x) c(paste0('coeff_',x),paste0('std_err_',x),paste0('t_',x),paste0('p_',x))))
    resultsTbl <- read.csv(text = paste(c('rsid',snp_names,cov_names,'RMSE'),collapse = ','))
    resultsTbl[1:length(rownames(dosageMatrix)),] <- NA
    for(k in 1:length(rownames(dosageMatrix))){
      resultsTbl[k,1] <- rownames(dosageMatrix)[k]
      resultsTbl[k,2:ncol(resultsTbl)] <- CalcInteractions(dosageMatrix[rownames(dosageMatrix)[k],],
                                                           dosageTarget,phenotypes,covariates)
    }
    data.table::fwrite(resultsTbl,file = paste0(path_out,'chunk',i,'.txt'),sep = '\t',col.names = T,row.names = F)
  },mc.cores = as.numeric(n_cores),ignore.interactive = T)
  #})
  saveRDS(allResultsTbl,file = paste0(path_out,'log.rds'))
}

##First read in the arguments listed at the command line
args=(commandArgs(TRUE))
print(args)
##args is now a list of character vectors
## First check to see if arguments are passed.
## Then cycle through each element of the list and evaluate the expressions.
if(length(args)==0){
  print("No arguments supplied.")
  #supply default values
  path <-  '~/bsu_scratch/LDL_Project_Data_Aug2019/Genotype_Data/'
  sample_file_prefix <- 'ukbb_LDL_metadata_with_PC'
  bgen_file_prefix <- 'ukb_imp_chr#_HRConly'
  chr <- '17'
  phenotype = 'LDLdirect'
  targetRS <- 'rs12916,rs17238484,rs5909,rs2303152,rs10066707,rs2006760'
  n_cores <- 1
  eur_only <- 1
  out_suffix <- 'test_full'
  cov <- 'sex,ages,bmi'
  PC <-  5
  MAF <- 0.05
  info <- 0.5
  med=1
  chunks <- '1,11'
  training_set <- 'NULL'
  if(out_suffix == ''){
    path_out <- paste0(gsub('Genotype_Data','Interaction_Data',path),gsub(',','_',targetRS),'_',phenotype,'/')
    path_out_chr <- paste0(gsub('Genotype_Data','Interaction_Data',path),gsub(',','_',targetRS),'_',phenotype,'/chr',chr,'/')
  }else{
    path_out <- paste0(gsub('Genotype_Data','Interaction_Data',path),gsub(',','_',targetRS),'_',phenotype,'_',out_suffix,'/')
    path_out_chr <- paste0(gsub('Genotype_Data','Interaction_Data',path),gsub(',','_',targetRS),'_',phenotype,'_',out_suffix,'/chr',chr,'/')
  }
  system(paste0('mkdir -p ',path_out))
  system(paste0('chmod a+rwx ',path_out))
  system(paste0('mkdir -p ',path_out_chr))
  system(paste0('chmod a+rwx ',path_out_chr))

}else if (length(args)!= 16){
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
       '# 15: Chunks\n',
       '# 16: Training Set\n',
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
           '# 14: Number of Cores:',
           '# 15: Chunks:',
           '# 16: Training Set:')
  variables <- c('path','sample_file_prefix','bgen_file_prefix','chr','phenotype','targetRS','out_suffix','eur_only',
                 'cov','PC','med','MAF','info','n_cores','chunks','training_set')
  for(i in 1:length(args)){
    eval(parse(text=paste0(variables[i],'=',"'",args[[i]],"'")))
    print(paste0(str[i],"   ",args[i]))
  }
  if(out_suffix == ''){
    path_out <- paste0(gsub('Genotype_Data','Interaction_Data',path),gsub(',','_',targetRS),'_',phenotype,'/')
    path_out_chr <- paste0(gsub('Genotype_Data','Interaction_Data',path),gsub(',','_',targetRS),'_',phenotype,'/chr',chr,'/')
  }else{
    path_out <- paste0(gsub('Genotype_Data','Interaction_Data',path),gsub(',','_',targetRS),'_',phenotype,'_',out_suffix,'/')
    path_out_chr <- paste0(gsub('Genotype_Data','Interaction_Data',path),gsub(',','_',targetRS),'_',phenotype,'_',out_suffix,'/chr',chr,'/')
  }
  system(paste0('mkdir -p ',path_out))
  system(paste0('chmod a+rwx ',path_out))
  system(paste0('mkdir -p ',path_out_chr))
  system(paste0('chmod a+rwx ',path_out_chr))
}
RunGxGInteractions(path,sample_file_prefix,bgen_file_prefix,chr,phenotype,targetRS,path_out_chr,eur_only,
                   cov,PC,med,MAF,info,n_cores,chunks,training_set)
