library(dplyr)
library(dbplyr)
library(data.table)
library(RSQLite)
interaction_path <- '~/bsu_scratch/LDL_Project_Data/Interaction_Data/HMGCR_LDL_known.txt'
interaction_results <- data.table::fread(interaction_path)
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
