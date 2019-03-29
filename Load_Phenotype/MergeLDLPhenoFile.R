library(data.table)
library(dplyr)
ldl_file <- data.table::fread('/home/zmx21/bsu_scratch/LDL_Project_Data/Genotype_Data/ldl_zmx2.csv')
pc_file <- data.table::fread('/home/zmx21/bsu_scratch/LDL_Project_Data/Genotype_Data/ukbb_metadata_with_PC.csv')

merged_df <- dplyr::left_join(dplyr::select(pc_file,c('UKB_genetic_ID',setdiff(colnames(pc_file),colnames(ldl_file))))
                              ,ldl_file,by = c("UKB_genetic_ID"="UKB_sample_ID"))

data.table::fwrite(merged_df,'/home/zmx21/bsu_scratch/LDL_Project_Data/Genotype_Data/ukbb_LDL_metadata_with_PC.csv')
