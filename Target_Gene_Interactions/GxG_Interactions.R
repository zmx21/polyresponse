RunGxGInteractions <- function(path,sample_file_prefix,bgen_file_prefix,phenotype,targetRS){
  chunkSize=100000
  
  library(dplyr)
  library(parallel)
  library(data.table)
  source('LoadBgen.R')
  
  #Find all rsids, and generate chunks to read.
  allRSIds <- FindAllRSIds(path,bgen_file_prefix) %>% dplyr::filter(rsid!=targetRS)
  rsIDChunks <- split(allRSIds,seq(nrow(allRSIds)-1)%/%chunkSize)
  samplesTbl <- LoadSamples(path,sample_file_prefix)
  
  #Intialize results vector
  resultsTbl <- data.frame('rsid'=rep('',nrow(allRSIds)),'p' = rep(NA,nrow(allRSIds)), 'coeff'=rep(NA,nrow(allRSIds)),stringsAsFactors = F)
  
  for(i in 1:length(rsIDChunks)){
    currentRSIdChunk <- rsIDChunks[[i]]
    genotype_data <- LoadBgen(path,bgen_file_prefix,currentRSIdChunk,targetRS)
    
    #Calculate allele dosage based on genotype probability
    alleleProbMatrix <- genotype_data$data
    remove(genotype_data)
    dosageMatrix <- matrix(0,nrow = nrow(alleleProbMatrix),ncol = ncol(alleleProbMatrix)) + alleleProbMatrix[,,'g=1'] + 2*alleleProbMatrix[,,'g=2']
    remove(alleleProbMatrix)
    
    #Construct Sample vs Phenotype Table
    samplePhenoTbl <- dplyr::select(samplesTbl,'samples'='ID_1',phenotype)
    colnames(dosageMatrix) <- samplePhenoTbl$samples
    
    #Function which calculates significance of interaction between two SNPs
    CalcInteractions <- function(dosageSubMatrix,phenotypes){
      #Merge into data frame
      data <- data.frame(SNP1=dosageSubMatrix[1,],SNP2=dosageSubMatrix[2,],Pheno=as.numeric(phenotypes))
      #Linear model, with GxG interaction plus individual main effect
      fit <- lm('Pheno~SNP1*SNP2',data = data)
      #Return significance of interaction term
      return(summary(fit)$coefficients[4,c(4,1)])
    }
    
    #save results
    nonTargetSnps <- setdiff(rownames(dosageMatrix),targetRS)
    for(k in 1:length(nonTargetSnps)){
      resultsTbl[((i-1)*chunkSize+k),1] <- nonTargetSnps[k]
      resultsTbl[((i-1)*chunkSize+k),c(2,3)] <- CalcInteractions(dosageMatrix[c(nonTargetSnps[k],targetRS),],samplePhenoTbl[,phenotype])
    }
  }
  data.table::fwrite(resultsTbl,file = paste0('~/',targetRS,'_',phenotype,'/',bgen_file_prefix,'.interactions.txt'),sep = '\t',row.names = F)
}
##First read in the arguments listed at the command line
args=(commandArgs(TRUE))

##args is now a list of character vectors
## First check to see if arguments are passed.
## Then cycle through each element of the list and evaluate the expressions.
if(length(args)==0){
  print("No arguments supplied.")
  ##supply default values
  path <- '~/blood_pressure_test/'
  sample_file_prefix <- 'ukbb_eur_all_sbp'
  bgen_file_prefix <- 'sbp_chr10'
  phenotype = 'sbp'
  targetRS <- 'rs603424'
  
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}
RunGxGInteractions(path,sample_file_prefix,bgen_file_prefix,phenotype,targetRS)
  