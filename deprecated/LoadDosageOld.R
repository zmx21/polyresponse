library(dplyr)
source('~/MRC_BSU_Internship/Load_Phenotype/Load_Phenotype.R')
source('~/MRC_BSU_Internship/Load_Bgen/LoadBgen.R')
source('~/MRC_BSU_Internship/SNP_Marginal_Effect/CalcMarginalEffect.R')
LoadDosage <- function(p_val_thresh,interaction_path,phenotype,r2_thresh,MAF,Info,targetRS){
  #Load rsid which passed interaction threshold.
  interaction_results <- data.table::fread(interaction_path) %>% dplyr::filter(p_int < p_val_thresh)
  
  #Connect to rsid annotation database. 
  anno_sql_name<- "all_snp_stats.sqlite"
  path <- '~/bsu_scratch/SQL/'
  curwd <- getwd()
  setwd(path)
  anno_con <- RSQLite::dbConnect(SQLite(), dbname = anno_sql_name)
  anno_db <- dplyr::tbl(anno_con,'all_snp_stats')
  #Get all rsids meeting criteria
  annotations <- anno_db %>% dplyr::filter(rsid %in% interaction_results$rsid) %>% collect()
  # annotations <- anno_db %>% dplyr::filter(rsid %in% interaction_results$rsid) %>% dplyr::select(rsid,minor_allele_frequency,info,chromosome) %>% collect()
  interaction_results <- interaction_results %>% dplyr::arrange(p_int) %>% dplyr::left_join(annotations,by=c('rsid'='rsid'))
  interaction_results <- interaction_results %>% dplyr::filter(as.numeric(minor_allele_frequency) > MAF & as.numeric(info) > Info) %>% dplyr::select(rsid,p_int,chromosome)
  interaction_results <- interaction_results %>% dplyr::distinct(rsid,.keep_all=T)
  #Load phenotype information and covariates
  path <-  '/mrc-bsu/scratch/zmx21/UKB_Data/'
  sample_file_prefix <- 'ukbb_metadata_with_PC'
  cov_names <- c('sex','ages','bmi')
  eur_only <- 1
  phenotypesAndCov <- LoadPhenotype(path,sample_file_prefix,phenotype,cov_names,eur_only=1,med=1)
  samplesToKeep <- phenotypesAndCov$samplesToKeep
  
  #Load dosage matrix of all inlcuded SNPS
  bgen_file_prefix <- 'ukb_imp_chr#_HRConly'
  dosageMatrix <- LoadBgen(path,bgen_file_prefix,interaction_results$rsid)
  dosageMatrix <- dosageMatrix[,samplesToKeep]
  
  #Generate independent set of predictors
  uniqueChr <- unique(interaction_results$chromosome)
  inclSet <- c()
  for(chr in uniqueChr){
    currentChr <- dplyr::filter(interaction_results,chromosome==chr)
    inclSet <- c(inclSet,currentChr$rsid[1])
    if(nrow(currentChr) > 1){
      LDMatrix <- cor(t(dosageMatrix[currentChr$rsid,])) ^ 2
      # LDMatrix <- GetLDMatrix(currentChr$rsid)
      currentChr <- dplyr::arrange(currentChr,p_int)
      for(i in 2:nrow(currentChr)){
        ldComparison = tryCatch({
          LDMatrix[currentChr$rsid[i],inclSet]
        }, warning = function(w) {
          ldComparison <- NA
        }, error = function(e) {
          ldComparison <- NA
        })
        if(any(is.na(ldComparison)) | any(ldComparison > r2_thresh)){
          next
        }
        inclSet <- c(inclSet,currentChr$rsid[i])
      }
    }
  }
  #only keep included set of SNPs in dosage matrix
  dosageMatrix <- dosageMatrix[inclSet,]
  
  #Name phenotypes and covariates
  names(phenotypesAndCov$phenotypes) <- colnames(dosageMatrix)
  rownames(phenotypesAndCov$covariates) <- colnames(dosageMatrix)
  
  
  #Dosage of target gene
  if(length(targetRS) < 2){
    dosageTarget <- LoadBgen(path,bgen_file_prefix,targetRS)
    dosageTarget <- dosageTarget[,samplesToKeep]
  }else{
    #Get marginal effects of rsid
    dosageVector <- LoadBgen(path,bgen_file_prefix,targetRS)
    dosageVector <- dosageVector[,samplesToKeep]
    
    beta_coeff <- CalcMarginalEffect(path,sample_file_prefix,bgen_file_prefix,phenotype,targetRS,
                                     eur_only,cov=paste(cov_names,collapse = ','),PC=5,med=1,1,F)$coeff
    #Flip alleles such that all are bp lowering
    dosageTarget <- dosageVector
    for(i in 1:nrow(dosageTarget)){
      if(beta_coeff[i] > 0){
        dosageTarget[i,] <- 2-dosageVector[i,]
      }
    }
    dosageTarget <- as.vector(abs(beta_coeff) %*% dosageTarget)
  }
  
  
  return(list(dosageTarget = dosageTarget,dosageMatrix = t(round(dosageMatrix,0)),
              phenotypes = phenotypesAndCov$phenotypes,covariates = phenotypesAndCov$covariates))
}
library(dplyr)
source('~/MRC_BSU_Internship/Load_Phenotype/Load_Phenotype.R')
source('~/MRC_BSU_Internship/Load_Bgen/LoadBgen.R')
source('~/MRC_BSU_Internship/SNP_Marginal_Effect/CalcMarginalEffect.R')
LoadDosage <- function(p_val_thresh,interaction_path,phenotype,r2_thresh,MAF,Info,targetRS){
  #Load rsid which passed interaction threshold.
  interaction_results <- data.table::fread(interaction_path) %>% dplyr::filter(p_int < p_val_thresh)
  
  #Connect to rsid annotation database. 
  anno_sql_name<- "all_snp_stats.sqlite"
  path <- '~/bsu_scratch/SQL/'
  curwd <- getwd()
  setwd(path)
  anno_con <- RSQLite::dbConnect(SQLite(), dbname = anno_sql_name)
  anno_db <- dplyr::tbl(anno_con,'all_snp_stats')
  #Get all rsids meeting criteria
  annotations <- anno_db %>% dplyr::filter(rsid %in% interaction_results$rsid) %>% collect()
  # annotations <- anno_db %>% dplyr::filter(rsid %in% interaction_results$rsid) %>% dplyr::select(rsid,minor_allele_frequency,info,chromosome) %>% collect()
  interaction_results <- interaction_results %>% dplyr::arrange(p_int) %>% dplyr::left_join(annotations,by=c('rsid'='rsid'))
  interaction_results <- interaction_results %>% dplyr::filter(as.numeric(minor_allele_frequency) > MAF & as.numeric(info) > Info) %>% dplyr::select(rsid,p_int,chromosome)
  interaction_results <- interaction_results %>% dplyr::distinct(rsid,.keep_all=T)
  #Load phenotype information and covariates
  path <-  '/mrc-bsu/scratch/zmx21/UKB_Data/'
  sample_file_prefix <- 'ukbb_metadata_with_PC'
  cov_names <- c('sex','ages','bmi')
  eur_only <- 1
  phenotypesAndCov <- LoadPhenotype(path,sample_file_prefix,phenotype,cov_names,eur_only=1,med=1)
  samplesToKeep <- phenotypesAndCov$samplesToKeep
  
  #Load dosage matrix of all inlcuded SNPS
  bgen_file_prefix <- 'ukb_imp_chr#_HRConly'
  dosageMatrix <- LoadBgen(path,bgen_file_prefix,interaction_results$rsid)
  dosageMatrix <- dosageMatrix[,samplesToKeep]
  
  #Generate independent set of predictors
  uniqueChr <- unique(interaction_results$chromosome)
  inclSet <- c()
  for(chr in uniqueChr){
    currentChr <- dplyr::filter(interaction_results,chromosome==chr)
    inclSet <- c(inclSet,currentChr$rsid[1])
    if(nrow(currentChr) > 1){
      LDMatrix <- cor(t(dosageMatrix[currentChr$rsid,])) ^ 2
      # LDMatrix <- GetLDMatrix(currentChr$rsid)
      currentChr <- dplyr::arrange(currentChr,p_int)
      for(i in 2:nrow(currentChr)){
        ldComparison = tryCatch({
          LDMatrix[currentChr$rsid[i],inclSet]
        }, warning = function(w) {
          ldComparison <- NA
        }, error = function(e) {
          ldComparison <- NA
        })
        if(any(is.na(ldComparison)) | any(ldComparison > r2_thresh)){
          next
        }
        inclSet <- c(inclSet,currentChr$rsid[i])
      }
    }
  }
  #only keep included set of SNPs in dosage matrix
  dosageMatrix <- dosageMatrix[inclSet,]
  
  #Name phenotypes and covariates
  names(phenotypesAndCov$phenotypes) <- colnames(dosageMatrix)
  rownames(phenotypesAndCov$covariates) <- colnames(dosageMatrix)
  
  
  #Dosage of target gene
  if(length(targetRS) < 2){
    dosageTarget <- LoadBgen(path,bgen_file_prefix,targetRS)
    dosageTarget <- dosageTarget[,samplesToKeep]
  }else{
    #Get marginal effects of rsid
    dosageVector <- LoadBgen(path,bgen_file_prefix,targetRS)
    dosageVector <- dosageVector[,samplesToKeep]
    
    beta_coeff <- CalcMarginalEffect(path,sample_file_prefix,bgen_file_prefix,phenotype,targetRS,
                                     eur_only,cov=paste(cov_names,collapse = ','),PC=5,med=1,1,F)$coeff
    #Flip alleles such that all are bp lowering
    dosageTarget <- dosageVector
    for(i in 1:nrow(dosageTarget)){
      if(beta_coeff[i] > 0){
        dosageTarget[i,] <- 2-dosageVector[i,]
      }
    }
    dosageTarget <- as.vector(abs(beta_coeff) %*% dosageTarget)
  }
  
  
  return(list(dosageTarget = dosageTarget,dosageMatrix = t(round(dosageMatrix,0)),
              phenotypes = phenotypesAndCov$phenotypes,covariates = phenotypesAndCov$covariates))
}
