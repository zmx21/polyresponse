LoadBgen <- function(path,bgen_file_prefix,rsIDs){
  library(rbgen)
  #Load bgen file
  curwd <- getwd()
  setwd(path)
  genotype_data <- rbgen::bgen.load(filename = paste0(bgen_file_prefix,'.bgen'),
                                    rsids = as.character(rsIDs$rsid))
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
  curwd <- getwd()
  setwd(path)
  
  bgenSummary <- system(paste0('bgenix -g ',bgen_file_prefix,'.bgen',' -list'),intern = T)
  #Remove comment headers
  bgenSummary <- bgenSummary[-c(1,length(bgenSummary))]
  bgenSummary <- data.table::fread(paste(bgenSummary,collapse = '\n'),sep = '\t',header = T,showProgress = F)
  rsIDs <- dplyr::select(bgenSummary,rsid,chromosome)
  setwd(curwd)
  
  return(rsIDs)
}
