library(dplyr)
int_results <- data.table::fread('allchr.txt')
snp_info <- data.table::fread('ukb_imp_allchr_HRConly_samplefilt_snpstats.txt')

int_results <- dplyr::select(int_results,rsid,coeff_int,std_err_int,t_int,p_int)
int_results <- dplyr::left_join(int_results,snp_info %>% dplyr::select(rsid,chromosome,position,alleleA,alleleB,alleleB_frequency,info),by=c('rsid'='rsid'))
int_results <- dplyr::arrange(int_results,chromosome,position)

#Excl HMGCR
chr <- 5
start_pos <- 74632154 - 2e6
end_pos <- 74657929 + 2e6
int_results <- int_results %>% dplyr::filter(!(chromosome == chr & position > start_pos & position < end_pos))

#Encode into required format
int_results <- dplyr::select(int_results,SNP = rsid,CHR = chromosome,POS = position, A1 = alleleA, A2 = alleleB, EAF = alleleB_frequency, Beta = coeff_int,se = std_err_int,P = p_int,INFO = info)
int_results <- dplyr::mutate(int_results,`A1/A2` = 'A1',N = 232419)
data.table::fwrite(int_results[,c('SNP','CHR','POS','A1','A2','A1/A2','EAF','Beta','se','P','INFO','N')],file = '~/Documents/Genetic_Epi_BSU_Project/gepi_22347_summary_stats.tsv',sep = '\t',quote = F,row.names = F,col.names = T)
