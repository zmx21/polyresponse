.libPaths(c("/home/zmx21/R/x86_64-redhat-linux-gnu-library/3.3"))
library(sqldf)
library(dplyr)
library(rbgen)
library(abind)
library(DBI)
library(RSQLite)
LoadSNPAnnotations <- function(sql_path,sql_db_name,rsid){
  curwd <- getwd()
  setwd(sql_path)
  con <- RSQLite::dbConnect(SQLite(), dbname = sql_db_name)
  #Get annotation info of rsid
  annotation <- do.call(rbind,lapply(rsid,function(x) RSQLite::dbGetQuery(con,paste0('SELECT * FROM all_snp_stats WHERE rsid=','"',x,'"'))))
  RSQLite::dbDisconnect(con)
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
    chr <- LoadSNPAnnotations(sql_path = "~/bsu_scratch/SQL/",sql_db_name = "all_snp_stats.sqlite",rsid = rsIDs) %>% dplyr::select(chromosome) %>% t() %>% as.numeric()
  }
  #Replace all wildcards with determined chromosomes
  bgen_file_prefix <- sapply(chr,function(x) gsub(pattern = '#',replacement = x,x = bgen_file_prefix))
  #Determine unique chr and load each chunk
  uniqueFiles <- unique(bgen_file_prefix)
  resultList <- vector(mode = 'list',length = length(uniqueFiles))
  
  for(i in 1:length(uniqueFiles)){
    currentFile <- uniqueFiles[i]
    currentrsIDs <- rsIDs[which(bgen_file_prefix == currentFile)]
    currentAlleleProbMatrix <- rbgen::bgen.load(filename = paste0(currentFile,'.bgen'),
                                      rsids = currentrsIDs)$data
    resultList[[i]] <- currentAlleleProbMatrix[,,'g=1'] + 2*currentAlleleProbMatrix[,,'g=2']
  }
  setwd(curwd)
  dosageMatrix <- do.call(rbind,resultList)
  return(dosageMatrix)
}
LoadSamples <- function(path,sample_file_prefix){
  #Load sample file
  samplesTbl <- read.table(file = paste0(path,sample_file_prefix,'.sample'),header = F,
                           stringsAsFactors = F)
  colnames(samplesTbl) <- samplesTbl[1,]#Add header as first row
  samplesTbl <- samplesTbl[-c(1,2),]#Remove first two rows
  return(samplesTbl)
}
# FindAllRSIds <- function(path,bgen_file_prefix){
#   library(dplyr)
#   library(data.table)
#   curwd <- getwd()
#   setwd(path)
#   
#   allFiles <- dir()
#   rsIDFile <- paste0(bgen_file_prefix,'_','rsid.txt')
#   #Run Bgenix tool if not previously run
#   if(!rsIDFile %in% allFiles){
#     bgenSummary <- system(paste0('bgenix -g ',bgen_file_prefix,'.bgen',' -list',' >',rsIDFile),intern = T)
#   }
#   bgenSummary <- data.table::fread(rsIDFile,sep = '\t',header = T,showProgress = F,skip = 1,fill = T)
#   bgenSummary <- bgenSummary[1:(nrow(bgenSummary)-1),]
#   rsIDs <- dplyr::select(bgenSummary,rsid,chromosome)
#   setwd(curwd)
#   
#   return(rsIDs)
# }

FindAllRSIds <- function(chr,MAF=-Inf,Info=-Inf){
  #Connect to rsid annotation database. 
  anno_sql_name<- "all_snp_stats.sqlite"
  path <- '~/bsu_scratch/SQL/'
  curwd <- getwd()
  setwd(path)
  anno_con <- RSQLite::dbConnect(SQLite(), dbname = anno_sql_name)
  anno_db <- dplyr::tbl(anno_con,'all_snp_stats')
  
  if(as.numeric(chr) < 10){
    chr <- paste0('0',chr)
  }
  
  #Get all rsids meeting criteria
  rsIDs <- anno_db %>% dplyr::filter(chromosome==chr & minor_allele_frequency > MAF & info > Info) %>% dplyr::select(rsid) %>% collect()
  setwd(curwd)
  RSQLite::dbDisconnect(anno_con)
  return(rsIDs)
}


# test <- LoadSNPAnnotations(sql_path = "~/bsu_scratch/SQL/",sql_db_name = "all_rsid.sqlite",rsid = 'rs193250821')
# test <- LoadBgen('~/bsu_scratch/UKB_Data/','ukb_imp_chr#_HRConly',c('rs193250821','rs537869321'))