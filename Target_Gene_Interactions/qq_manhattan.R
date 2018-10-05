library(qqman)
library(RSQLite)
library(dplyr)
library(pbmcapply)
#Connect to rsid annotation database. 
anno_sql_name<- "all_snp_stats.sqlite"
path <- '~/bsu_scratch/SQL/'
setwd(path)
anno_con <- RSQLite::dbConnect(SQLite(), dbname = anno_sql_name)
anno_db <- tbl(anno_con,'all_snp_stats')


#Read in results
results_tbl <- data.table::fread(paste0('~/bsu_scratch/UKB_Data/rs1262894_sbp_eur_only/chr10.txt'),header = F)
colnames(results_tbl) <- c('rsid','p','coeff','rsq')
results_tbl$p_bonf <- p.adjust(results_tbl$p_int,method = 'bonferroni')

sig_results <- results_tbl %>% dplyr::filter(p_bonf<0.2)
#Join with annotation information
sig_results_annotation <- anno_db %>% dplyr::select(rsid,chromosome,position) %>% dplyr::filter(rsid %in% sig_results$rsid) %>% collect()
sig_results <- dplyr::left_join(sig_results,sig_results_annotation,by=c('rsid'='rsid'))

#Mahattan plot
manhattan(data.frame(BP=as.numeric(sig_results$position),CHR=as.numeric(sig_results$chromosome),P=sig_results$p_bonf,SNP=sig_results$rsid),suggestiveline = F,genomewideline = -log10(0.05),annotateTop = T,ylab='-log((',annotatePval = 0.0001)


#QQPlot
#Convert p values to chi-sq.
PToChiSq <- function(p,df){
  return(qt(p/2,df) ^2)
}
# results_tbl_anno <- anno_db %>% dplyr::select(rsid,chromosome) %>% dplyr::filter(rsid %in% results_tbl$rsid) %>% collect()

library(snpStats)
results_tbl$chi_sq <- unlist(pbmclapply(results_tbl$p,function(x) PToChiSq(x,460209),mc.cores = 30),use.names = F)
qq.chisq(results_tbl$chi_sq,main = 'rs1262894 chr10 EUR (Age,Sex,BMI)')

maf_info_filtered <-dplyr::filter(anno_db,rsid %in% results_tbl$rsid & minor_allele_frequency > 0.05 & info > 0.5) %>% dplyr::select(rsid) %>% collect() %>% dplyr::left_join(results_tbl,by=c('rsid'='rsid'))
qq.chisq(maf_info_filtered$chi_sq,main = 'rs1262894 chr10 EUR MAF and Info Filtered')

