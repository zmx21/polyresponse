library(dplyr)
chr <- 10
print(chr)
path <- '~/bsu_scratch/UKB_Data/rs603424_sbp/'
interactionResults <- data.table::fread(paste0(path,'chr',chr,'.txt'))
colnames(interactionResults) <- c('rsid','p','coeff')

SNPdb <- data.table::fread(paste0('~/bsu_scratch/All_VCF/ensembl/homo_sapiens-chr',chr,'.vcf'))
interactionResults <- dplyr::left_join(interactionResults,SNPdb,by=c('rsid'='ID'))
