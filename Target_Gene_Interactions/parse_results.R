####################################################################################
#Parse GxG interaction results, keeping only rsid and p value.
####################################################################################
library(dplyr)
library(data.table)
parse_results <- function(input,output){
  all_snps <- data.table::fread(input,header = T,select=c('rsid','p_int'))
  #Remove missing rows
  parsed_snps <- all_snps %>% dplyr::filter(!is.na(p_int))
  parsed_snps <- parsed_snps %>% dplyr::distinct(rsid,.keep_all=T)
  #write results
  data.table::fwrite(parsed_snps,file=output,row.names = F,sep='\t')
}
