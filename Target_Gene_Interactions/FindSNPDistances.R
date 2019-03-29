library(RSQLite)
library(dplyr)
interaction_results <- data.table::fread('~/parsed_interaction/CACNA1D_sbp.txt') 
filt_int_results <- interaction_results %>% dplyr::filter(p_int<4e-5)
#Connect to rsid annotation database. 
anno_sql_name<- "all_snp_stats.sqlite"
path <- '~/bsu_scratch/SQL/'
curwd <- getwd()
setwd(path)
anno_con <- RSQLite::dbConnect(SQLite(), dbname = anno_sql_name)
anno_db <- dplyr::tbl(anno_con,'all_snp_stats')
filt_int_results <-  dplyr::filter(anno_db,rsid %in% filt_int_results$rsid) %>% collect() %>% dplyr::left_join(filt_int_results,by=c('rsid'='rsid'))
saveRDS(filt_int_results,file='~/CACNA1D_4e-5.rds')

known_snps <- c('rs4652519','rs1949872','rs12494691','rs588076','rs597446','rs541458','rs7196683','rs2429427','rs1029765','rs2246709','rs2740574','rs776746','rs1123617','rs3745009','rs1137617','rs12946454','rs11191548','rs12346562','rs1104514')
known_interactions <- dplyr::filter(anno_db,rsid %in% known_snps) %>% collect() %>%  dplyr::left_join(interaction_results,by=c('rsid'='rsid'))
saveRDS(known_interactions,file='~/CACNA1D_Known.rds')