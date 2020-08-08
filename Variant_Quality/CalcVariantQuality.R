####################################################################################
#Calcuates quality metrics (Info Score, MAF etc ) of variants in a bgen file,
#using BASH script called "calc_variant_quality", located in "path" argument.
#The BASH script calls QCTOOL.
#
#Input: chromosome to calculate metrics, assuming bgen files are organized by chr.
#Output: files written in "./variant" within "path", containing output from qctool.
####################################################################################
source('../Load_Bgen/LoadBgen.R')
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
  allRSIds <- FindAllRSIds(chr)
  allRSIds <- unique(allRSIds$rsid)
  chunkSize <- ceiling(length(allRSIds) / n_cores)
  rsIDChunks <- split(allRSIds,seq(length(allRSIds)-1)%/%chunkSize)
  #Write file which contains the rsid chunks (used by QC tools)
  for(i in 1:length(rsIDChunks)){
    write(rsIDChunks[[i]],file = paste0(path_out,'chunk',i,'_rsid.txt'))
  }
  #Call "calc_variant_quality", a BASH script which calls QCTOOL
  mclapply(1:length(rsIDChunks),function(i) system(paste0('./calc_variant_quality ',chr,' ',paste0('./variants/chr',chr,'/','chunk',i,'_rsid.txt'),' ',paste0('./variants/chr',chr,'/','chunk',i,'_stats.txt'))),mc.cores = n_cores)
}
path <-  '/mrc-bsu/scratch/zmx21/UKB_Data/' #Path of data
n_cores <- 16
args=(commandArgs(TRUE))
chr <- args[1] #Chr to process
CalcVariantQuality(path,chr,n_cores)