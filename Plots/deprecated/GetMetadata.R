
library(RcppEigen)
library(dplyr)
library(parallel)
library(data.table)
library(pbmcapply)

source('~/MRC_BSU_Internship_LDL/Load_Phenotype/Load_Phenotype.R')
source('~/MRC_BSU_Internship_LDL/Load_Bgen/LoadBgen.R')
source('~/MRC_BSU_Internship_LDL/SNP_Phenotype_Association/CalcSnpPhenoAssociation.R')




path <-  '~/bsu_scratch/LDL_Project_Data/Genotype_Data/'
sample_file_prefix <- 'ukbb_LDL_metadata_with_PC'
bgen_file_prefix <- 'ukb_imp_chr#_HRConly'
chr <- '5'
targetRS <- 'rs12916,rs17238484,rs5909,rs2303152,rs10066707,rs2006760'
phenotype <- 'LDLdirect'
cov <- 'sex,ages,bmi,sbp,dbp'
PC <-  5
med <- 1
eur_only <- 1
n_cores <- 1

#Parse covariate argument
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
phenotypesAndCov <- LoadPhenotype(path,sample_file_prefix,phenotype,cov_names,eur_only,med,F,c('sbp','dbp'))
phenotypes <- phenotypesAndCov$phenotypes
covariates <- as.matrix(phenotypesAndCov$covariates)
samplesToKeep <- phenotypesAndCov$samplesToKeep


#Split rsid based on comma
targetRS <- unlist(strsplit(targetRS,split = ','))
dosageTarget <- LoadBgen(path,bgen_file_prefix,targetRS)

#Calculate gene score if more than 1 target SNP listed, otherwise load single SNP dosage vector.
if(length(targetRS) < 2){
  #Load single SNP dosage vector
  dosageTarget <- LoadBgen(path,bgen_file_prefix,targetRS)
  #keep samples without missing information
  dosageTarget <- dosageTarget[,samplesToKeep]

}else{
  #Get dosage vectors of all SNPs
  dosageVector <- LoadBgen(path,bgen_file_prefix,targetRS)
  #keep samples without missing information
  dosageVector <- dosageVector[,samplesToKeep]

  #Calculate main effects of each individual SNP
  beta_coeff <- CalcSnpPhenoAssociation(path,sample_file_prefix,bgen_file_prefix,phenotype,targetRS,as.numeric(eur_only),cov=cov,PC=as.numeric(PC),med=as.numeric(med),as.numeric(n_cores),F)$coeff
  #Flip alleles such that all are LDL lowering
  dosageTarget <- dosageVector
  for(i in 1:nrow(dosageTarget)){
    if(beta_coeff[i] > 0){
      dosageTarget[i,] <- 2-dosageVector[i,]
    }
  }
  dosageTarget <- as.vector(abs(beta_coeff) %*% dosageTarget)
}
save.image('~/bsu_scratch/LDL_Project_Data/Metadata.RData')
CalcStats <- function(covariates,dosageTarget,phenotypes,lipdbin){
  if(!all.equal(nrow(covariates),length(dosageTarget),length(phenotypes))){
    return(NULL)
  }

  n_samples <- length(dosageTarget)
  mean_age <- mean(covariates[,'ages'])
  sd_age <- sd(covariates[,'ages'])
  mean_bmi <- mean(covariates[,'bmi'])
  sd_bmi <- sd(covariates[,'bmi'])
  num_female <- length(covariates[,'sex']) - sum(covariates[,'sex'])
  perc_female <- num_female / n_samples
  mean_sbp <- mean(covariates[,'sbp'],na.rm = T)
  sd_sbp <- sd(covariates[,'sbp'],na.rm = T)
  mean_dbp <- mean(covariates[,'dbp'],na.rm = T)
  sd_dbp <- sd(covariates[,'dbp'],na.rm = T)

  mean_LDL <- mean(phenotypes)
  sd_LDL <- sd(phenotypes)


  perc_med <- sum(lipdbin == 'Current') / length(lipdbin)

  mean_score <- mean(dosageTarget)
  sd_score <- sd(dosageTarget)
  return(data.frame(n_samples=n_samples,
                    age=mean_age,
                    sd_age=sd_age,
                    mean_bmi=mean_bmi,
                    sd_bmi=sd_bmi,
                    perc_female=perc_female,
                    mean_LDL=mean_LDL,
                    sd_LDL=sd_LDL,
                    mean_score=mean_score,
                    sd_score=sd_score,
                    mean_sbp=mean_sbp,
                    sd_sbp=sd_sbp,
                    mean_dbp=mean_dbp,
                    sd_dbp=sd_dbp,
                    perc_med=perc_med))
}

#Calculate stats for overall data
overall_stats <- CalcStats(covariates = covariates,dosageTarget = dosageTarget,phenotypes = phenotypes,lipdbin = phenotypesAndCov$lipdbin)
below_med <- which(dosageTarget <= median(dosageTarget))
above_med <- which(dosageTarget > median(dosageTarget))
below_med_stats <- CalcStats(covariates = covariates[below_med,],dosageTarget = dosageTarget[below_med],phenotypes = phenotypes[below_med],lipdbin = phenotypesAndCov$lipdbin[below_med])
above_med_stats <- CalcStats(covariates = covariates[above_med,],dosageTarget = dosageTarget[above_med],phenotypes = phenotypes[above_med],lipdbin = phenotypesAndCov$lipdbin[above_med])


all_stats <- t(rbind(overall_stats,below_med_stats,above_med_stats))
colnames(all_stats) <- c('Overall','HMGCR Score <= Median','HMGCR Score > Median')
all_stats <- signif(all_stats,4)

p_age <- t.test(covariates[below_med,'ages'],covariates[above_med,'ages'])$p.value
sex_cont_tbl <- t(table(c(rep('Below Med',length(below_med)),rep('Above Med',length(above_med))),c(covariates[below_med,'sex'],covariates[above_med,'sex'])))
p_sex <- chisq.test(sex_cont_tbl)$p.value
p_bmi <- t.test(covariates[below_med,'bmi'],covariates[above_med,'bmi'])$p.value
med_cont_tbl <- t(table(c(rep('Below Med',length(below_med)),rep('Above Med',length(above_med))),c(phenotypesAndCov$lipdbin[below_med],phenotypesAndCov$lipdbin[above_med])))
p_med <- chisq.test(med_cont_tbl)$p.value
p_sbp <- t.test(covariates[below_med,'sbp'],covariates[above_med,'sbp'])$p.value
p_dbp <- t.test(covariates[below_med,'dbp'],covariates[above_med,'dbp'])$p.value
p_LDL <- t.test(phenotypes[below_med],phenotypes[above_med])$p.value
p_score <- t.test(dosageTarget[below_med],dosageTarget[above_med])$p.value

p_val_df <- c(n_samples=NA,p_age=p_age,p_bmi=p_bmi,p_sex=p_sex,p_LDL=p_LDL,p_score=p_score,p_sbp=p_sbp,p_dbp=p_dbp,p_med=p_med)
p_val_df <- as.data.frame(signif(p_val_df,4))
colnames(p_val_df) <- 'p-value'

