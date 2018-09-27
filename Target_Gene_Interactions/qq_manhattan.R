library(qqman)
library(RSQLite)
library(dplyr)
library(pbmcapply)
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
colnames(results_tbl) <- c('rsid','p','coeff')
results_tbl$p_bonf <- p.adjust(results_tbl$p,method = 'bonferroni')

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
# results_tbl$chi_sq <- unlist(pbmclapply(results_tbl$p,function(x) PToChiSq(x,460209),mc.cores = 30),use.names = F)

library(snpStats)
results_sample <- results_tbl[sample(1:nrow(results_tbl),1000000),]
results_sample$chi_sq <- unlist(pbmclapply(results_sample$p,function(x) PToChiSq(x,460209),mc.cores = 30),use.names = F)

qq.chisq(results_sample$chi_sq)
