library(sqldf)
library(dplyr)
library(rbgen)
library(abind)
LoadSNPAnnotations <- function(sql_path,sql_db_name,rsid){
  curwd <- getwd()
  setwd(sql_path)
  con <- RSQLite::dbConnect(SQLite(), dbname = sql_db_name)
  #Get annotation info of rsid
  annotation <- do.call(rbind,lapply(rsid,function(x) RSQLite::dbGetQuery(con,paste0('SELECT * FROM all_rsid WHERE rsid=','"',x,'"'))))
  #Remove duplicated rsid (different mutation at same pos)
  annotation <- dplyr::distinct(annotation,position,chromosome,.keep_all=T)
  setwd(curwd)
  return(annotation)
}

#Import should be directory where bgen files are, and the prefix of the files (# represents chr number)
LoadBgen <- function(path,bgen_file_prefix,rsIDs,chr=NULL){
  #Load bgen file
  curwd <- getwd()
  setwd(path)
  #If chr information of each snp is not provided, look up annotation to determine chromosome.
  if(is.null(chr)){
    chr <- LoadSNPAnnotations(sql_path = "~/bsu_scratch/SQL/",sql_db_name = "all_rsid.sqlite",rsid = rsIDs) %>% dplyr::select(chromosome) %>% t() %>% as.numeric()
  }
  #Replace all wildcards with determined chromosomes
  bgen_file_prefix <- sapply(chr,function(x) gsub(pattern = '#',replacement = x,x = bgen_file_prefix))
  #Determine unique chr and load each chunk
  uniqueFiles <- unique(bgen_file_prefix)
  resultList <- vector(mode = 'list',length = length(uniqueFiles))
  
  for(i in 1:length(uniqueFiles)){
    currentFile <- uniqueFiles[i]
    currentrsIDs <- rsIDs[which(bgen_file_prefix == currentFile)]
    resultList[[i]] <- rbgen::bgen.load(filename = paste0(currentFile,'.bgen'),
                                      rsids = currentrsIDs)$data
  }
  setwd(curwd)
  alleleProbMatrix <- abind(resultList,along = 1)
  return(alleleProbMatrix)
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

# test <- LoadSNPAnnotations(sql_path = "~/bsu_scratch/SQL/",sql_db_name = "all_rsid.sqlite",rsid = 'rs193250821')
# test <- LoadBgen('~/bsu_scratch/UKB_Data/','ukb_imp_chr#_HRConly',c('rs193250821','rs537869321'))