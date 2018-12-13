library(dplyr)
source('~/MRC_BSU_Internship/Gene_Phenotype_Association/GetSNPsOfGene.R')
source('~/MRC_BSU_Internship/Gene_Phenotype_Association/GenePhenotypeAssociation.R')
CACNA1D_SNPs <- CollectBetaCoeff('CACNA1D',10000,10000,'sbp',0.01,0.5,16)  %>% dplyr::filter(!is.na('P.value'))
includedSNPS <- c('rs3821843','rs7340705','rs113210396','rs312487','rs11719824','rs3774530','rs3821856')
CACNA1D_SNPs <- CACNA1D_SNPs %>% dplyr::filter(rsid %in% includedSNPS)

#Connect to annotation database
anno_sql_name<- "all_snp_stats.sqlite"
path <- '~/bsu_scratch/SQL/'
setwd(path)
anno_con <- RSQLite::dbConnect(SQLite(), dbname = anno_sql_name)
anno_db <- tbl(anno_con,'all_snp_stats')
anno_db <- anno_db %>% dplyr::filter(rsid %in% includedSNPS) %>% collect()

CACNA1D_SNPs <- dplyr::left_join(anno_db,CACNA1D_SNPs,by=c('rsid' = 'rsid')) %>% dplyr::select(rsid,'P-value'=p,'Main Effect'=coeff,'Chr'=chromosome,'Position'=position,'Minor Allele'=minor_allele,'Major Allele'=major_allele,'MAF' = minor_allele_frequency)
CACNA1D_SNPs$`Main Effect` <- signif(CACNA1D_SNPs$`Main Effect`,4)
CACNA1D_SNPs$`P-value` <- signif(CACNA1D_SNPs$`P-value`,4)
CACNA1D_SNPs$MAF <- signif(as.numeric(CACNA1D_SNPs$MAF),4)

data.table::fwrite(CACNA1D_SNPs[order(CACNA1D_SNPs$`P-value`),],file='~/figures/final_figures/CACNA1D_SNP_stats.csv',sep = ',')
