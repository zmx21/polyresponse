library(pbmcapply)
CollectBetaCoeff <- function(gene_name,upstream_dist,downstream_dist,phenotype,MAF,info,n_cores){
  source('~/MRC_BSU_Internship/Target_Score/get_snps_of_gene.R')
  source('~/MRC_BSU_Internship/SNP_Marginal_Effect/CalcMarginalEffect.R')
  #Obtain all SNPs which map to the gene
  print('Finding SNPs')
  all_snps <- as.data.frame(AllSNPsOfGene(gene_name,upstream_dist,downstream_dist))
  all_snps$minor_allele_frequency <- as.numeric(all_snps$minor_allele_frequency)
  all_snps$info <- as.numeric(all_snps$info)
  all_snps <- all_snps %>% dplyr::filter(minor_allele_frequency > !!MAF & info > !!info)
  
  #Calculate marginal effect of all snps using the UKBB data
  path <-  '~/bsu_scratch/UKB_Data/'
  sample_file_prefix <- 'ukbb_metadata_with_PC'
  bgen_file_prefix <- 'ukb_imp_chr#_HRConly'

  UKBB_marginal_effect <- CalcMarginalEffect(path,sample_file_prefix,bgen_file_prefix,phenotype,all_snps$rsid,1,cov='sex,ages,bmi',PC=5,med=1,n_cores,F)
  
  return(UKBB_marginal_effect)
}
# PKD2L1_SNPs <- CollectBetaCoeff('PKD2L1',1000,500,'sbp',0.05,0.5,10)

# SCNN1A_SNPs <- CollectBetaCoeff('SCNN1A',1000,500,10)
# SCNN1B_SNPs <- CollectBetaCoeff('SCNN1B',1000,500,10)
# SCNN1D_SNPs <- CollectBetaCoeff('SCNN1D',100000,100000,0.05,0.5,30)
# SCNN1G_SNPs <- CollectBetaCoeff('SCNN1G',1000,500,10)
