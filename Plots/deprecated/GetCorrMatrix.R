source('~/MRC_BSU_Internship_LDL/Load_Bgen/LoadBgen.R')
includedSNPS <- c('rs12916','rs17238484','rs5909','rs2303152','rs10066707','rs2006760')
# interactionResults <- data.table::fread('~/bsu_scratch/LDL_Project_Data/HMGCR_Interaction_Stats.csv')
# includedSNPS <- interactionResults[sapply(interactionResults$genes,function(x) grepl('IFNA',x)),]$rsid

path <-  '~/bsu_scratch/LDL_Project_Data/Genotype_Data/'
bgen_file_prefix <- 'ukb_imp_chr#_HRConly'

dosageMatrix <- LoadBgen(path,bgen_file_prefix,includedSNPS)
r2_matrix <- matrix(nrow = length(includedSNPS),ncol = length(includedSNPS))
row.names(r2_matrix) <- includedSNPS
colnames(r2_matrix) <- includedSNPS

for(i in 1:nrow(r2_matrix)){
  for(j in i:ncol(r2_matrix)){
    r2_matrix[i,j] <- cor(as.vector(dosageMatrix[i,]),as.vector(dosageMatrix[j,])) ^ 2
  }
}

library(corrplot)
corrplot(r2_matrix,type='upper',method = 'color',addCoef.col = 'black')
