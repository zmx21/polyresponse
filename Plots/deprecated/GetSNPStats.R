library(dplyr)
source('~/MRC_BSU_Internship_LDL/SNP_Phenotype_Association/CalcSnpPhenoAssociation.R')
path <-  '~/bsu_scratch/LDL_Project_Data/Genotype_Data/'
sample_file_prefix <- 'ukbb_LDL_metadata_with_PC'
bgen_file_prefix <- 'ukb_imp_chr#_HRConly'
phenotype <- 'LDLdirect'
n_cores <- 1

includedSNPS <- c('rs12916','rs17238484','rs5909','rs2303152','rs10066707','rs2006760')
HMGCR_SNPs <- CalcSnpPhenoAssociation(path,sample_file_prefix,bgen_file_prefix,phenotype,includedSNPS,1,cov='sex,ages,bmi',PC=5,med=1,n_cores,F)

#Connect to annotation database
anno_sql_name<- "all_snp_stats.sqlite"
path <- '~/bsu_scratch/SQL/'
setwd(path)
anno_con <- RSQLite::dbConnect(SQLite(), dbname = anno_sql_name)
anno_db <- tbl(anno_con,'all_snp_stats')
anno_db <- anno_db %>% dplyr::filter(rsid %in% includedSNPS) %>% collect()

HMGCR_SNPs <- dplyr::left_join(anno_db,HMGCR_SNPs,by=c('rsid' = 'rsid')) %>% dplyr::select(rsid,'p-value'=p,'Beta'=coeff,'Chr'=chromosome,'Position'=position,'Minor Allele'=minor_allele,'Major Allele'=major_allele,'MAF' = minor_allele_frequency)
HMGCR_SNPs$`Beta` <- signif(HMGCR_SNPs$`Beta`,4)
HMGCR_SNPs$`p-value` <- signif(HMGCR_SNPs$`p-value`,4)
HMGCR_SNPs$MAF <- signif(as.numeric(HMGCR_SNPs$MAF),4)

data.table::fwrite(HMGCR_SNPs[order(HMGCR_SNPs$`p-value`),],file='~/bsu_scratch/LDL_Project_Data/HMGCR_SNP_know_stats.csv',sep = ',')


