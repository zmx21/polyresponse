CalcVariantQuality <- function(path,chr,n_cores){
  path_out <- paste0(path,'variants/','chr',chr,'/')
  system(paste0('mkdir -p ',path_out))
  setwd(path)
  library(dplyr)
  library(parallel)
  library(data.table)
  library(pbmcapply)
  #Find all rsids, and generate chunks to read.
  print('Loading rsID')
  source('~/MRC_BSU_Internship/Load_Bgen/LoadBgen.R')
  allRSIds <- FindAllRSIds(chr)
  allRSIds <- unique(allRSIds$rsid)
  chunkSize <- ceiling(length(allRSIds) / n_cores)
  rsIDChunks <- split(allRSIds,seq(length(allRSIds)-1)%/%chunkSize)
  
  for(i in 1:length(rsIDChunks)){
    write(rsIDChunks[[i]],file = paste0(path_out,'chunk',i,'_rsid.txt'))
  }
  
  mclapply(1:length(rsIDChunks),function(i) system(paste0('./calc_variant_quality ',chr,' ',paste0('./variants/chr',chr,'/','chunk',i,'_rsid.txt'),' ',paste0('./variants/chr',chr,'/','chunk',i,'_stats.txt'))),mc.cores = n_cores)
}
path <-  '/mrc-bsu/scratch/zmx21/UKB_Data/'
n_cores <- 16
args=(commandArgs(TRUE))
chr <- args[1]
CalcVariantQuality(path,chr,n_cores)