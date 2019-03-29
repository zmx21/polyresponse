####################################################################################
#Identifies a set of independent SNPs within given gene region
#Input:Target gene name, upstream/downstream distance, phenotype, p-value threshold, and r2 threshold
#Output: List of independent SNPs, with their invididual marignal effect
####################################################################################
source('~/MRC_BSU_Internship_LDL/Load_Bgen/LoadBgen.R')
source('~/MRC_BSU_Internship_LDL/Load_Phenotype/Load_Phenotype.R')
source('~/MRC_BSU_Internship_LDL/Gene_Phenotype_Association/GetSNPsOfGene.R')
source('~/MRC_BSU_Internship_LDL/SNP_Phenotype_Association/CalcSnpPhenoAssociation.R')

library(RcppEigen)
library(pbmcapply)
library(dplyr)

#Collects beta coefficients of all SNPs within a gene region, and calculates mean effect of each individual SNP
CollectBetaCoeff <- function(gene_name,upstream_dist,downstream_dist,phenotype,MAF,info,n_cores){
  #Obtain all SNPs which map to the gene
  print('Finding SNPs')
  all_snps <- as.data.frame(AllSNPsOfGene(gene_name,upstream_dist,downstream_dist))
  all_snps$minor_allele_frequency <- as.numeric(all_snps$minor_allele_frequency)
  all_snps$info <- as.numeric(all_snps$info)
  all_snps <- all_snps %>% dplyr::filter(minor_allele_frequency > !!MAF & info > !!info)
  
  #Calculate mean effect of all snps using the UKBB data
  path <-  '~/bsu_scratch/LDL_Project_Data/Genotype_Data/'
  sample_file_prefix <- 'ukbb_LDL_metadata_with_PC'
  bgen_file_prefix <- 'ukb_imp_chr#_HRConly'
  
  UKBB_mean_effects <- CalcSnpPhenoAssociation(path,sample_file_prefix,bgen_file_prefix,phenotype,all_snps$rsid,1,cov='sex,ages,bmi',PC=5,med=1,n_cores,F)
  
  return(UKBB_mean_effects)
}

#Generate a set of indenpendent SNPs within a gene region, ranked according to associate with phenotype.
IterativePruning <- function(gene_name,phenotype,upstream_dist,downstream_dist,p_thresh,r2_thresh,n_cores){
  betaCoeff <- CollectBetaCoeff(gene_name,upstream_dist,downstream_dist,phenotype,0.01,0.5,n_cores)
  betaCoeffFilt <- dplyr::filter(betaCoeff,p < p_thresh)
  if(nrow(betaCoeffFilt) < 2){
    print('Less than 2 sig SNPs')
    if(nrow(betaCoeffFilt) == 0){
      betaCoeffFilt <- dplyr::arrange(betaCoeff,p) %>% {.[1,]}
    }
    return(list(rsid=betaCoeffFilt$rsid,coeff = betaCoeffFilt))
    
  }
  print('Finding LD')
  #Load phenotype information and covariates
  path <-  '~/bsu_scratch/LDL_Project_Data/Genotype_Data/'
  sample_file_prefix <- 'ukbb_LDL_metadata_with_PC'
  cov_names <- c('sex','ages','bmi','PC1','PC2','PC3','PC4','PC5')
  eur_only <- 1
  phenotypesAndCov <- LoadPhenotype(path,sample_file_prefix,phenotype,cov_names,eur_only=1,med=1)
  samplesToKeep <- phenotypesAndCov$samplesToKeep
  
  #Load dosage matrix of all inlcuded SNPS
  bgen_file_prefix <- 'ukb_imp_chr#_HRConly'
  dosageMatrix <- LoadBgen(path,bgen_file_prefix,betaCoeffFilt$rsid)
  dosageMatrix <- dosageMatrix[,samplesToKeep]
  
  LDMatrix <- cor(t(dosageMatrix)) ^ 2
  betaCoeffFilt <- dplyr::arrange(betaCoeffFilt,p)
  inclSet <- betaCoeffFilt$rsid[1]
  for(i in 2:nrow(betaCoeffFilt)){
    ldComparison = tryCatch({
      LDMatrix[betaCoeffFilt$rsid[i],inclSet]
    }, warning = function(w) {
      ldComparison <- NA
    }, error = function(e) {
      ldComparison <- NA
    })
    if(any(is.na(ldComparison)) | any(ldComparison > r2_thresh)){
      next
    }
    inclSet <- c(inclSet,betaCoeffFilt$rsid[i])
  }
  return(list(rsid=inclSet,coeff = dplyr::filter(betaCoeffFilt,rsid%in%!!inclSet)))
}

#Run iterative pruning.
target_names <- c('HMGCR')
target_LDL <- lapply(target_names,function(x){
  result <- IterativePruning(x,'LDLdirect',10000,10000,5e-8,0.3,16)
  cbind(data.frame(gene_name=rep(x,nrow(result$coeff)),stringsAsFactors = F),result$coeff)})
target_LDL <- dplyr::left_join(do.call(rbind,target_LDL),AllSNPsOfGene('HMGCR',10000,10000),by = c('rsid'='rsid'))
data.table::fwrite(target_LDL,row.names = F,file = '~/bsu_scratch/LDL_Project_Data/Genotype_Data/HMGCR_SNPs.csv')

#Get effect sizes for known snps
known_snps = c('rs12916','rs17238484','rs5909','rs2303152','rs10066707','rs2006760')
rs12916_Fit <- CalcSnpPhenoAssociation(path = '~/bsu_scratch/LDL_Project_Data/Genotype_Data/',
                                            sample_file_prefix = 'ukbb_LDL_metadata_with_PC',
                                            bgen_file_prefix = 'ukb_imp_chr#_HRConly',
                                            phenotype = 'LDLdirect',
                                            known_snps[1],1,cov='sex,ages,bmi',PC=5,med=1,1,T)
saveRDS(rs12916_Fit,'~/bsu_scratch/LDL_Project_Data/Genotype_Data/rs12916_fit.rds')


target_LDL_known <- CalcSnpPhenoAssociation(path = '~/bsu_scratch/LDL_Project_Data/Genotype_Data/',
                                            sample_file_prefix = 'ukbb_LDL_metadata_with_PC',
                                            bgen_file_prefix = 'ukb_imp_chr#_HRConly',
                                            phenotype = 'LDLdirect',
                                            known_snps,1,cov='sex,ages,bmi',PC=5,med=1,16,F)
target_LDL_known <- dplyr::left_join(target_LDL_known,AllSNPsOfGene('HMGCR',100000,100000),by = c('rsid'='rsid'))
data.table::fwrite(target_LDL_known,row.names = F,file = '~/bsu_scratch/LDL_Project_Data/Genotype_Data/HMGCR_SNPs_Known.csv')



####################################################################################
#Calcuate the association between gene score and phenotype
#Input: vector of rsid and their respective beta coefficients.
#Output: fit result of the gene score and its association with phenotype 
####################################################################################


#Called by CalcMeanEffectGeneScore, constructs a linear model with/without covariate
#to measure genotype-phenotype beta parameter and significance. 
LinearFitGeneScore <- function(dosageSubMatrix,phenotypes,covariates){
  #Create model matrix, with main and interaction effects of two SNPs. 
  if(ncol(covariates) > 0){
    mdlMat <- cbind('Intercept' = rep(1,length(dosageSubMatrix)),'gene_score' = dosageSubMatrix,covariates)
  }else{
    mdlMat <- cbind('Intercept' = rep(1,length(dosageSubMatrix)),'gene_score' = dosageSubMatrix)
  }
  
  fit <- RcppEigen::fastLm(y = phenotypes,X = mdlMat)
  #Return coefficient and significance of interaction term
  results <- c(summary(fit)$coefficients[2,c(1,4)],summary(fit)$r.squared)
  names(results) <- c('coeff','p','rsq')
  return(results)
}

#Calculates mean effect of gene score with phenotype
CalcMeanEffectGeneScore <- function(path,sample_file_prefix,bgen_file_prefix,phenotype,rsid,eur_only,cov,beta_coeff,med=T,PC=5){
  #Decide what columns to load based on what covariates were specificed
  if(cov!=''){
    cov_names <- unlist(strsplit(x=cov,split = ','))
  }else{
    cov_names <- c()
  }
  #Add PC to covariates, if specified
  if(PC > 0){
    cov_names <- c(cov_names,sapply(1:PC,function(x) paste0('PC',x)))
  }
  
  #Load Phenotype and Covariates
  print('Loading Phenotypes')
  phenotypesAndCov <- LoadPhenotype(path,sample_file_prefix,phenotype,cov_names,eur_only,med)
  phenotypes <- phenotypesAndCov$phenotypes
  covariates <- phenotypesAndCov$covariates
  samplesToKeep <- phenotypesAndCov$samplesToKeep
  
  dosageVector <- LoadBgen(path,bgen_file_prefix,rsid)
  dosageVector <- dosageVector[,samplesToKeep]
  
  dosageVector <- beta_coeff %*% dosageVector
  fit <- LinearFitGeneScore(as.vector(dosageVector),phenotypes,covariates)
  return(fit)
}
