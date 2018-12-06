####################################################################################
#Loads genotype dosage of specified SNPs, in the format which would be used to build RF
#Input: P-value threshold and path to where interactions results are already stored
#MAF, Info, and r2 threshold should also be provided. 
#Output: dataframe containing genotype of target and mediating SNPs, phenotype, and covariates
####################################################################################

library(dplyr)
source('~/MRC_BSU_Internship/Load_Phenotype/Load_Phenotype.R')
source('~/MRC_BSU_Internship/Load_Bgen/LoadBgen.R')
source('~/MRC_BSU_Internship/SNP_Phenotype_Association/CalcSnpPhenoAssociation.R')
source('~/MRC_BSU_Internship/Target_Gene_Interactions/CalcInteractions.R')
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
  dbDisconnect(anno_con)
  interaction_results <- interaction_results %>% dplyr::arrange(p_int) %>% dplyr::left_join(annotations,by=c('rsid'='rsid'))
  interaction_results <- interaction_results %>% dplyr::filter(as.numeric(minor_allele_frequency) > MAF & as.numeric(info) > Info) %>% dplyr::select(rsid,p_int,chromosome)
  interaction_results <- interaction_results %>% dplyr::distinct(rsid,.keep_all=T)
  #Load phenotype information and covariates
  path <-  '/mrc-bsu/scratch/zmx21/UKB_Data/'
  sample_file_prefix <- 'ukbb_metadata_with_PC'
  cov_names <- c('sex','ages','bmi')
  eur_only <- 1
  
  #Load phenotype and samples to keep
  phenotypesAndCov <- LoadPhenotype(path,sample_file_prefix,phenotype,c(cov_names,'PC1','PC2','PC3','PC4','PC5'),eur_only=1,med=1)
  samplesToKeep <- phenotypesAndCov$samplesToKeep
  
  #Dosage of target gene
  bgen_file_prefix <- 'ukb_imp_chr#_HRConly'
  if(length(targetRS) < 2){
    dosageTarget <- LoadBgen(path,bgen_file_prefix,targetRS)
    dosageTarget <- dosageTarget[,samplesToKeep]
  }else{
    #Get mean effects of rsid
    dosageVector <- LoadBgen(path,bgen_file_prefix,targetRS)
    dosageVector <- dosageVector[,samplesToKeep]
    
    beta_coeff <- CalcSnpPhenoAssociation(path,sample_file_prefix,bgen_file_prefix,phenotype,targetRS,
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
  
  #Generate independent set of predictors
  uniqueChr <- unique(interaction_results$chromosome)
  inclSet <- c()
  inclSetDosage <- matrix(nrow = sum(samplesToKeep),ncol = 0)
  pb = txtProgressBar()
  for(chr in uniqueChr){
    # setTxtProgressBar(pb,value = which(chr == uniqueChr)/length(uniqueChr))
    currentChr <- dplyr::filter(interaction_results,chromosome==chr) %>% dplyr::arrange(p_int)
    inclSet <- c(inclSet,currentChr$rsid[1])
    currentChrDosage <- as.matrix(LoadBgen(path,bgen_file_prefix,currentChr$rsid[1])[,samplesToKeep])
    if(min(dim(currentChrDosage)) > 1){
      names(dosageTarget) <- rownames(currentChrDosage)
      interactionSig <- apply(currentChrDosage,1,function(x) CalcInteractions(as.vector(x),dosageTarget,phenotypesAndCov$phenotypes,as.matrix(phenotypesAndCov$covariates))[12])
      currentChrDosage <- as.matrix(currentChrDosage[which.min(interactionSig),])
    }
    if(nrow(currentChr) > 1){
      for(i in 2:nrow(currentChr)){
        currentRS <- currentChr$rsid[i]
        currentDosage <- LoadBgen(path,bgen_file_prefix,currentChr$rsid[i])
        currentDosage <- as.matrix(currentDosage[,samplesToKeep])
        if(min(dim(currentDosage)) > 1){
          names(dosageTarget) <- rownames(currentChrDosage)
          interactionSig <- apply(currentDosage,1,function(x) CalcInteractions(as.vector(x),dosageTarget,phenotypesAndCov$phenotypes,as.matrix(phenotypesAndCov$covariates))[12])
          currentDosage <- as.matrix(currentDosage[which.min(interactionSig),])
        }
        #Calculate correlation with every snp in the included set
        if(ncol(currentChrDosage)==1){
          pairwiseLD <- cor(as.vector(t(currentDosage)),as.vector(t(currentChrDosage))) ^ 2
        }else{
          pairwiseLD <- cor(currentDosage,currentChrDosage) ^ 2
        }
        #don't include snp if exceed correlation threshold
        if(any(pairwiseLD > r2_thresh)){
          next
        }
        inclSet <- c(inclSet,currentRS)
        currentChrDosage <- cbind(currentChrDosage,currentDosage)
      }
    }
    inclSetDosage <- cbind(inclSetDosage,currentChrDosage)
  }
  dosageMatrix <- inclSetDosage
  colnames(dosageMatrix) <- inclSet
  #Name phenotypes and covariates
  names(phenotypesAndCov$phenotypes) <- rownames(dosageMatrix)
  rownames(phenotypesAndCov$covariates) <- rownames(dosageMatrix)
  names(dosageTarget) <- rownames(dosageMatrix)
  return(list(dosageTarget = dosageTarget,dosageMatrix = round(dosageMatrix,0),
              phenotypes = phenotypesAndCov$phenotypes,covariates = phenotypesAndCov$covariates))
}