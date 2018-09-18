LoadBgen <- function(path,bgen_file_prefix,rsIDs){
  library(rbgen)
  #Load bgen file
  curwd <- getwd()
  setwd(path)
  genotype_data <- rbgen::bgen.load(filename = paste0(bgen_file_prefix,'.bgen'),
                                    rsids = rsIDs)
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
  library(data.table)
  curwd <- getwd()
  setwd(path)
  
  allFiles <- dir()
  rsIDFile <- paste0(bgen_file_prefix,'_','rsid.txt')
  #Run Bgenix tool if not previously run
  if(!rsIDFile %in% allFiles){
    bgenSummary <- system(paste0('bgenix -g ',bgen_file_prefix,'.bgen',' -list',' >',rsIDFile),intern = T)
  }
  bgenSummary <- data.table::fread(rsIDFile,sep = '\t',header = T,showProgress = F,skip = 1,fill = T)
  bgenSummary <- bgenSummary[1:(nrow(bgenSummary)-1),]
  rsIDs <- dplyr::select(bgenSummary,rsid,chromosome)
  setwd(curwd)
  
  return(rsIDs)
}
