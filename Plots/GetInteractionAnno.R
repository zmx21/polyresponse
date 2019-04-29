library(dplyr)
library(dbplyr)
library(data.table)
library(RSQLite)
interaction_path <- '~/bsu_scratch/LDL_Project_Data/Interaction_Data/HMGCR_LDL_known.txt'
#interaction_path <- '~/bsu_scratch/LDL_Project_Data/Interaction_Data/PCSK9_LDL_known.txt'

interaction_results <- data.table::fread(interaction_path)

WriteAnnovarInput <- function(){
  #Filter variants based on threshold
  interaction_results <- dplyr::filter(interaction_results,p_int < 1e-2) %>% dplyr::arrange(p_int)
  
  #Connect to rsid annotation database. 
  anno_sql_name<- "all_snp_stats.sqlite"
  path <- '~/bsu_scratch/SQL/'
  curwd <- getwd()
  setwd(path)
  anno_con <- RSQLite::dbConnect(SQLite(), dbname = anno_sql_name)
  anno_db <- dplyr::tbl(anno_con,'all_snp_stats')
  #Get all rsids meeting criteria
  anno <- anno_db %>% dplyr::filter(rsid %in% interaction_results$rsid) %>% collect()
  dbDisconnect(anno_con)
  interaction_results <- interaction_results %>% dplyr::left_join(anno,by=c('rsid'='rsid')) 
  chr_pos <- interaction_results %>% dplyr::select(chromosome,position)
  chr_pos <- cbind(chr_pos,data.frame(position2=chr_pos$position))
  output <- cbind(chr_pos,dplyr::select(interaction_results,alleleA,alleleB,p_int))
  output$p_int <- as.character(output$p_int)
  
  data.table::fwrite(output,'~/bsu_scratch/LDL_Project_Data/Interaction_Data/annovar/HMGCR_int_parsed',col.names = F,sep='\t')
  #data.table::fwrite(output,'~/bsu_scratch/LDL_Project_Data/Interaction_Data/annovar/PCSK9_int_parsed',col.names = F,sep='\t')
}

GetAllInteractionMetadata <- function(){
  #Filter variants based on threshold
  interaction_results <- dplyr::filter(interaction_results,p_int < 5e-5) %>% dplyr::arrange(p_int)
  
  #Connect to rsid annotation database. 
  anno_sql_name<- "all_snp_stats.sqlite"
  path <- '~/bsu_scratch/SQL/'
  curwd <- getwd()
  setwd(path)
  anno_con <- RSQLite::dbConnect(SQLite(), dbname = anno_sql_name)
  anno_db <- dplyr::tbl(anno_con,'all_snp_stats')
  
  #Get all rsids meeting criteria
  anno <- anno_db %>% dplyr::filter(rsid %in% interaction_results$rsid) %>% collect()
  dbDisconnect(anno_con)
  interaction_results <- interaction_results %>% dplyr::left_join(anno,by=c('rsid'='rsid')) 
  interaction_results <- interaction_results %>% dplyr::select(rsid,'p-value'=p_int,'Chr'=chromosome,'Position'=position,'Minor Allele'=minor_allele,'Major Allele'=major_allele,'MAF' = minor_allele_frequency)
  interaction_results$`p-value` <- signif(interaction_results$`p-value`,4)
  interaction_results$MAF <- signif(as.numeric(interaction_results$MAF),4)
  

  #Read in annovar results
  resultPath <- '~/bsu_scratch/LDL_Project_Data/Interaction_Data/annovar/'
  file_name <- 'HMGCR_int_parsed.variant_function'
  annovarResult <- data.table::fread(paste0(resultPath,file_name))
  annovarResultFilt <- annovarResult %>% dplyr::filter(V8 < 5e-5) %>% dplyr::arrange(V8) 
  annovarResultFilt <- annovarResultFilt %>% dplyr::select(type=V1,genes=V2,chr=V3,pos=V4,p=V8)

  #Seperate according to comma
  allGenes <- annovarResultFilt$genes
  #Remove none entries and NM entries
  removeEntries <- sapply(allGenes,function(x) grepl('NONE',x)) | sapply(allGenes,function(x) grepl('NM_',x))
  allGenes[removeEntries] <- sapply(allGenes[removeEntries],function(x) gsub(x = x,pattern = 'NONE(dist=NONE),',replacement = '',fixed = T))
  allGenes[removeEntries] <- sapply(allGenes[removeEntries],function(x) gsub(x = x,pattern = 'NONE(dist=NONE)',replacement = '',fixed = T))
  allGenes[removeEntries] <- sapply(allGenes[removeEntries],function(x) unlist(strsplit(x,'(',fixed = T))[1],USE.NAMES = F)
  allGenes <- sapply(allGenes,function(x) gsub(x=x,pattern = ',',replacement = '\\'))
  
  interaction_results$Gene <- allGenes
  
  data.table::fwrite(interaction_results[order(interaction_results$`p-value`),],file='~/bsu_scratch/LDL_Project_Data/HMGCR_Interaction_Stats.csv',sep = ',')
  
  
}














