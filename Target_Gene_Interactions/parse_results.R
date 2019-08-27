
####################################################################################
#Parse GxG interaction results, keeping only rsid and p value.
####################################################################################
source('~/MRC_BSU_Internship_LDL/Gene_Phenotype_Association/GetSNPsOfGene.R')
library(dplyr)
library(data.table)
library(qqman)

parse_results <- function(input,output,gene_region,verbose = F){
  all_snps <- data.table::fread(input,header = T,select=c('rsid','p_int'))
  #Remove missing rows
  parsed_snps <- all_snps %>% dplyr::filter(!is.na(p_int))
  dup_snps <- parsed_snps$rsid[duplicated(parsed_snps$rsid)]
  parsed_snps <- parsed_snps %>% dplyr::distinct(rsid,.keep_all=T)
  
  HMGCR_snps <- AllSNPsOfGene('HMGCR',2e6,2e6) %>% 
    dplyr::left_join(parsed_snps,by=c('rsid'='rsid')) %>% dplyr::filter(!is.na(p_int))
  if(verbose == T){
    HMGCR_broad_snps <- AllSNPsOfGene('HMGCR',5e6,5e6) %>% 
      dplyr::left_join(parsed_snps,by=c('rsid'='rsid')) %>% dplyr::filter(!is.na(p_int))
    ggplot(data.frame(pos=as.numeric(HMGCR_broad_snps$position),P=HMGCR_broad_snps$p_int)) 
    + aes(x=pos,y=-log10(P)) + geom_point()
  }
  #write results
  parsed_snps <- dplyr::filter(parsed_snps,!rsid %in% HMGCR_snps$rsid)
  data.table::fwrite(parsed_snps,file=paste0(output,'.txt'),row.names = F,sep='\t')
  write(dup_snps,sep = '\n',file = paste0(output,'.dup','.txt'))
}
