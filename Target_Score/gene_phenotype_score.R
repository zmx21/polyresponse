library(pbmcapply)
CollectBetaCoeff <- function(gene_name,upstream_dist,downstream_dist,n_cores){
  source('~/MRC_BSU_Internship/Target_Score/get_snps_of_gene.R')
  source('~/MRC_BSU_Internship/SNP_Marginal_Effect/CalcMarginalEffect.R')
  #Obtain all SNPs which map to the gene
  all_snps <- AllSNPsOfGene(gene_name,upstream_dist,downstream_dist)
  all_snps <- all_snps %>% dplyr::filter(grepl("rs",rsid))
  
  #Calculate GWAS P values
  GWAS_Results <- data.table::fread('~/bsu_scratch/GWAS/Ehret_2016_summstats.txt')
  GWAS_Results <- dplyr::inner_join(all_snps,GWAS_Results,by=c('rsid'='rsid'))
  
  #Calculate marginal effect of all snps using the UKBB data
  path <-  '/mrc-bsu/scratch/zmx21/UKB_Data/'
  sample_file_prefix <- 'ukbb_eur_all_sbp'
  bgen_file_prefix <- 'ukb_imp_chr#_HRConly'
  phenotype <- 'sbp'
  samplesTbl <- LoadSamples(path,sample_file_prefix)
  
  UKBB_marginal_effect <- unlist(pbmclapply(GWAS_Results$rsid,function(x) CalcMarginalEffect(path,samplesTbl,bgen_file_prefix,phenotype,x),mc.cores = n_cores))
  
  
  return(cbind(GWAS_Results,data.frame(UKBB_marginal_effect=UKBB_marginal_effect)))
}
# PKD2L1_SNPs <- CollectBetaCoeff('PKD2L1',1000,500,10)

# SCNN1A_SNPs <- CollectBetaCoeff('SCNN1A',1000,500,10)
# SCNN1B_SNPs <- CollectBetaCoeff('SCNN1B',1000,500,10)
SCNN1D_SNPs <- CollectBetaCoeff('SCNN1D',1000,500,10)
# SCNN1G_SNPs <- CollectBetaCoeff('SCNN1G',1000,500,10)
