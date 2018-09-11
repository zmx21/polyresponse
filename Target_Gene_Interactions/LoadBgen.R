LoadBgen <- function(path,bgen_file_prefix,rsIDs,targetRS){
  library(rbgen)
  #Load bgen file
  curwd <- getwd()
  setwd(path)
  genotype_data <- rbgen::bgen.load(filename = paste0(bgen_file_prefix,'.bgen'),
                                    rsids = c(as.character(rsIDs$rsid),targetRS))
  setwd(curwd)
  return(genotype_data)
}
LoadSamples <- function(path,sample_file_prefix){
  #Load sample file
  samplesTbl <- read.table(file = paste0(path,sample_file_prefix,'.sample'),header = F,
                           stringsAsFactors = F)
  colnames(samplesTbl) <- samplesTbl[1,]#Add header as first row
  samplesTbl <- samplesTbl[-c(1,2),]#Remove first two rows
  return(samplesTbl)
}
FindAllRSIds <- function(path,bgen_file_prefix){
  library(dplyr)
  bgenSummary <- system(paste0('bgenix -g ',path,bgen_file_prefix,'.bgen',' -list'),intern = T)
  #Remove comment headers
  bgenSummary <- bgenSummary[-c(1,length(bgenSummary))]
  bgenSummary <- read.delim(text = bgenSummary,header = T)
  rsIDs <- dplyr::select(bgenSummary,rsid,chromosome)
  return(rsIDs)
}