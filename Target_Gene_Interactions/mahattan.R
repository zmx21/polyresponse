library(qqman)
library(RSQLite)
library(dplyr)

#Connect to rsid annotation database. 
anno_sql_name<- "all_rsid.sqlite"
path <- '~/bsu_scratch/SQL/'
setwd(path)
anno_con <- RSQLite::dbConnect(SQLite(), dbname = anno_sql_name)
anno_db <- tbl(anno_con,'all_rsid')


#Read in results
target <- "rs603424"
phenotype <- 'sbp'
results_tbl <- data.table::fread(paste0('~/bsu_scratch/UKB_Data/',target,'_',phenotype,'/all_chr.txt'))
results_tbl$p_bonf <- p.adjust(results_tbl$p,method = 'bonferroni')
results_tbl$p_fdr <- p.adjust(results_tbl$p,method = 'fdr')

sig_results <- results_tbl %>% dplyr::filter(p_fdr<0.2)
#Join with annotation information
sig_results_annotation <- anno_db %>% dplyr::select(rsid,chromosome,position) %>% dplyr::filter(rsid %in% sig_results$rsid) %>% collect()
sig_results <- dplyr::left_join(sig_results,sig_results_annotation,by=c('rsid'='rsid'))
saveRDS(sig_results,'~/bsu_scratch/UKB_Data/rs603424_sbp/sig_results.rds')

manhattan(data.frame(BP=as.numeric(sig_results$position),CHR=as.numeric(sig_results$chromosome),P=sig_results$p_fdr,SNP=sig_results$rsid),suggestiveline = F,genomewideline = -log10(0.05),annotateTop = T,ylab='-log(adj_p)',annotatePval = 0.0001)
