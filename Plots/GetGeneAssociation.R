library(dplyr)
source('~/MRC_BSU_Internship_LDL/SNP_Phenotype_Association/CalcSnpPhenoAssociation.R')
source('~/MRC_BSU_Internship_LDL/Gene_Phenotype_Association/GenePhenotypeAssociation.R')

path <-  '~/bsu_scratch/LDL_Project_Data/Genotype_Data/'
sample_file_prefix <- 'ukbb_LDL_metadata_with_PC'
bgen_file_prefix <- 'ukb_imp_chr#_HRConly'
phenotype <- 'LDLdirect'
n_cores <- 1

includedSNPS <- c('rs12916','rs17238484','rs5909','rs2303152','rs10066707','rs2006760')
HMGCR_SNPs <- CalcSnpPhenoAssociation(path,sample_file_prefix,bgen_file_prefix,phenotype,includedSNPS,1,cov='sex,ages,bmi',PC=5,med=1,n_cores,F)
beta_coeff <- as.numeric(HMGCR_SNPs$coeff)

mean_effect_fit <- CalcMeanEffectGeneScore(path = path,
                                           sample_file_prefix = sample_file_prefix,
                                           bgen_file_prefix = bgen_file_prefix,
                                           phenotype = phenotype,
                                           rsid = includedSNPS,
                                           eur_only = 1,
                                           cov = 'ages,sex,bmi',
                                           beta_coeff = beta_coeff,
                                           med=T,
                                           PC=5,
                                           verbose=T)