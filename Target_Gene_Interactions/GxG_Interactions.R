library(RcppEigen)
#Function which calculates significance of interaction between two SNPs
CalcInteractions <- function(dosageSubMatrix,dosageTarget,phenotypes){
  #Create model matrix, with interaction between SNPs.
  mdlMat <- cbind('Intercept' = rep(1,length(dosageSubMatrix)),'InteractingSNP' = dosageSubMatrix,'TargetSNP' = dosageTarget,'Interaction' = dosageSubMatrix * dosageTarget)
  fit <- RcppEigen::fastLm(y = phenotypes,X = mdlMat)
  # #Merge into data frame
  # data <- data.frame(SNP1=dosageSubMatrix,SNP2=dosageTarget,Pheno=as.numeric(phenotypes))
  # #Linear model, with GxG interaction plus individual main effect
  # fit <- lm('Pheno~SNP1*SNP2',data = data)
  #Return significance of interaction term
  return(summary(fit)$coefficients[4,c(4,1)])
}

RunGxGInteractions <- function(path,sample_file_prefix,bgen_file_prefix,chr,phenotype,targetRS,target_chr,path_out,n_cores){
  chunkSize=200
  
  library(dplyr)
  library(parallel)
  library(data.table)
  library(pbmcapply)
  source('LoadBgen.R')
  
  #Genotype info of target gene
  target_file_prefix <- gsub(pattern = '#',replacement = target_chr,x = bgen_file_prefix)
  dosageTarget <- LoadBgen(path,target_file_prefix,data.frame(rsid=targetRS))
  dosageTarget <- dosageTarget$data
  dosageTarget <- matrix(0,nrow = nrow(dosageTarget),ncol = ncol(dosageTarget)) + dosageTarget[,,'g=1'] + 2*dosageTarget[,,'g=2']
  
  #Generate prefix
  bgen_file_prefix <- gsub(pattern = '#',replacement = chr,x = bgen_file_prefix)
  
  #Find all rsids, and generate chunks to read.
  print('Loading rsID')
  allRSIds <- FindAllRSIds(path,bgen_file_prefix) %>% dplyr::filter(rsid!=targetRS) %>% dplyr::distinct(rsid)
  rsIDChunks <- split(allRSIds,seq(nrow(allRSIds)-1)%/%chunkSize)
  print('Loading Samples')
  samplesTbl <- LoadSamples(path,sample_file_prefix)
  
  #Run chunks in parallel
  print('Calculating Interactions')
  allResultsTbl <- pbmclapply(rsIDChunks,function(currentRSIdChunk) {
    genotype_data <- LoadBgen(path,bgen_file_prefix,currentRSIdChunk)
    
    #Calculate allele dosage based on genotype probability
    alleleProbMatrix <- genotype_data$data
    remove(genotype_data)
    dosageMatrix <- matrix(0,nrow = nrow(alleleProbMatrix),ncol = ncol(alleleProbMatrix)) + alleleProbMatrix[,,'g=1'] + 2*alleleProbMatrix[,,'g=2']
    remove(alleleProbMatrix)
    
    #Construct Sample vs Phenotype Table
    samplePhenoTbl <- dplyr::select(samplesTbl,'samples'='ID_1',phenotype)
    phenotypes<- as.numeric(samplePhenoTbl[,phenotype])
    #save results
    nonTargetSnps <- setdiff(rownames(dosageMatrix),targetRS)
    dosageMatrix <- dosageMatrix[,!is.nan(phenotypes)]
    dosageTarget <- dosageTarget[,!is.nan(phenotypes)]
    phenotypes <- phenotypes[!is.nan(phenotypes)]
    
    #Intialize results vector
    resultsTbl <- data.frame('rsid'=rep('',length(nonTargetSnps)),'p' = rep(NA,length(nonTargetSnps)), 'coeff'=rep(NA,length(nonTargetSnps)),stringsAsFactors = F)
    for(k in 1:length(nonTargetSnps)){
      resultsTbl[k,1] <- nonTargetSnps[k]
      resultsTbl[k,c(2,3)] <- CalcInteractions(dosageMatrix[nonTargetSnps[k],],dosageTarget,phenotypes)
    }
    return(resultsTbl)
  },mc.cores = as.numeric(n_cores),ignore.interactive = T)
  data.table::fwrite(do.call(rbind,allResultsTbl),file = paste0(path_out,bgen_file_prefix,'.interactions.txt'),sep = '\t',row.names = F)
}
##First read in the arguments listed at the command line
args=(commandArgs(TRUE))
##args is now a list of character vectors
## First check to see if arguments are passed.
## Then cycle through each element of the list and evaluate the expressions.
if(length(args)==0){
  print("No arguments supplied.")
  #supply default values
  # path <-  '~/rds/hpc-work/UKB_Data/'
  # sample_file_prefix <- 'ukbb_eur_all_sbp'
  # bgen_file_prefix <- 'ukb_imp_chr#_HRConly'
  # chr <- '10'
  # phenotype = 'sbp'
  # targetRS <- 'rs603424'
  # target_chr <- '10'
  # path_out <- paste0(path,targetRS,'_',phenotype)
  # n_cores=2
  path <-  '/home/zmx21/blood_pressure_test/'
  sample_file_prefix <- 'ukbb_eur_all_sbp'
  bgen_file_prefix <- 'sbp_chr#'
  chr <- '10'
  phenotype = 'sbp'
  targetRS <- 'rs603424'
  target_chr <- '10'
  path_out <- paste0(path,targetRS,'_',phenotype)
  n_cores=2
  system(paste0('mkdir -p ',path_out))
}else if (length(args)!=8){
  stop("You need to supply:\n",
       "# 1: Input Path\n",
       '# 2: Sample File Prefix\n',
       "# 3: Bgen File Rrefix\n",
       "# 4: Chromosome to Analyze\n",
       "# 5: Phenotype Identifier (column name in .sample file) \n",
       "# 6: rsid of Target SNP\n",
       "# 7: Chromosome of The Target SNP\n",
       '# 8: Number of Cores\n',
       "Exiting...", call.=FALSE)
  
}else{
  print('All Arguments Supplied')
  str <- c("# 1: Input Path:",
           '# 2: Sample File Prefix:',
           "# 3: Bgen File Prefix:",
           "# 4: Chromosome to Analyze:",
           "# 5: Phenotype Identifier (column name in .sample file):",
           "# 6: rsid of Target SNP:",
           "# 7: Chromosome of The Target SNP:",
           '# 8: Number of Cores:')
  variables <- c('path','sample_file_prefix','bgen_file_prefix','chr','phenotype','targetRS','target_chr','n_cores')
  for(i in 1:length(args)){
    eval(parse(text=paste0(variables[i],'=',"'",args[[i]],"'")))
    print(paste0(str[i],"   ",args[i]))
  }
  path_out <- paste0(path,targetRS,'_',phenotype,'/')
  system(paste0('mkdir -p ',path_out))
}
RunGxGInteractions(path,sample_file_prefix,bgen_file_prefix,chr,phenotype,targetRS,target_chr,path_out,n_cores)
