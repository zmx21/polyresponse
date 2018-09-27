library(RcppEigen)
#Function which calculates significance of interaction between two SNPs.
CalcInteractions <- function(dosageSubMatrix,dosageTarget,phenotypes){
  #Create model matrix, with main and interaction effects of two SNPs. 
  mdlMat <- cbind('Intercept' = rep(1,length(dosageSubMatrix)),'InteractingSNP' = dosageSubMatrix,'TargetSNP' = dosageTarget,'Interaction' = dosageSubMatrix * dosageTarget)
  fit <- RcppEigen::fastLm(y = phenotypes,X = mdlMat)
  #Return coefficient and significance of interaction term
  return(c(summary(fit)$coefficients[4,c(4,1)],summary(fit)$r.squared))
}

RunGxGInteractions <- function(path,sample_file_prefix,bgen_file_prefix,chr,phenotype,targetRS,target_chr,path_out,n_cores){
  chunkSize=50
  
  library(dplyr)
  library(parallel)
  library(data.table)
  library(pbmcapply)
  source('/home/zmx21/MRC_BSU_Internship/Load_Bgen/LoadBgen.R')
  
  #Genotype info of target gene
  dosageTarget <- LoadBgen(path,bgen_file_prefix,targetRS,target_chr)

  #Find all rsids, and generate chunks to read.
  print('Loading rsID')
  allRSIds <- FindAllRSIds(path,gsub(pattern = '#',replacement = chr,x = bgen_file_prefix)) %>% dplyr::filter(rsid!=targetRS)
  allRSIds <- unique(allRSIds$rsid)
  rsIDChunks <- split(allRSIds,seq(length(allRSIds)-1)%/%chunkSize)
  print('Loading Samples')
  samplesTbl <- LoadSamples(path,sample_file_prefix)
  
  #Run chunks in parallel
  print(paste0('Calculating Interactions: ',length(rsIDChunks),' chunks'))
  allResultsTbl <- pbmclapply(1:length(rsIDChunks),function(i) {
    currentRSIdChunk <- rsIDChunks[[i]]
    dosageMatrix <- LoadBgen(path,bgen_file_prefix,currentRSIdChunk,rep(chr,length(currentRSIdChunk)))
    
    #Construct Sample vs Phenotype Table
    samplePhenoTbl <- dplyr::select(samplesTbl,'samples'='ID_1',phenotype)
    phenotypes<- as.numeric(samplePhenoTbl[,phenotype])
    #save results
    nonTargetSnps <- setdiff(rownames(dosageMatrix),targetRS)
    dosageMatrix <- dosageMatrix[,!is.nan(phenotypes)]
    dosageTarget <- dosageTarget[,!is.nan(phenotypes)]
    phenotypes <- phenotypes[!is.nan(phenotypes)]
    
    #Intialize results vector
    resultsTbl <- data.frame('rsid'=rep('',length(nonTargetSnps)),'p' = rep(NA,length(nonTargetSnps)), 'coeff'=rep(NA,length(nonTargetSnps)),'rsq'=rep(NA,length(nonTargetSnps)),stringsAsFactors = F)
    for(k in 1:length(nonTargetSnps)){
      resultsTbl[k,1] <- nonTargetSnps[k]
      resultsTbl[k,c(2,3,4)] <- CalcInteractions(dosageMatrix[nonTargetSnps[k],],dosageTarget,phenotypes)
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
  path <-  '/mrc-bsu/scratch/zmx21/UKB_Data/blood_pressure_test/'
  sample_file_prefix <- 'ukbb_eur_all_sbp'
  bgen_file_prefix <- 'sbp_chr10'
  chr <- '10'
  phenotype = 'sbp'
  targetRS <- 'rs603424'
  target_chr <- '10'
  n_cores=10
  path_out <- paste0(path,targetRS,'_',phenotype,'/')
  path_out_chr <- paste0(path,targetRS,'_',phenotype,'/chr',chr,'/')
  system(paste0('mkdir -p ',path_out))
  system(paste0('chmod a+rwx ',path_out))
  system(paste0('mkdir -p ',path_out_chr))
  system(paste0('chmod a+rwx ',path_out_chr))

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
  path_out_chr <- paste0(path,targetRS,'_',phenotype,'/chr',chr,'/')
  system(paste0('mkdir -p ',path_out))
  system(paste0('chmod a+rwx ',path_out))
  system(paste0('mkdir -p ',path_out_chr))
  system(paste0('chmod a+rwx ',path_out_chr))
}
RunGxGInteractions(path,sample_file_prefix,bgen_file_prefix,chr,phenotype,targetRS,target_chr,path_out_chr,n_cores)
system(paste0('cat ', path_out_chr,'* > ',path_out_chr,'chr',chr,'.txt'))
