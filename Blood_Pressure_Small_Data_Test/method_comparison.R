library(dplyr)
library(data.table)
library(rbgen)
#Obtain reults from PLINK and manually calculated using R
plinkResults <- data.table::fread('/mrc-bsu/scratch/zmx21/UKB_Data/blood_pressure_test/plink.epi.qt')
manualResults <- data.table::fread('/mrc-bsu/scratch/zmx21/UKB_Data/blood_pressure_test/rs603424_sbp/chr10/chr10.txt')

#Sort both results according to rsID
plinkResults <- plinkResults %>% dplyr::select(rsid=SNP2,p_plink=P,coeff_plink=BETA_INT) 
manualResults <- manualResults %>% dplyr::select(rsid=V1,p_manual=V2,coeff_manual=V3) 

#Load genotype data to look at imputation information
curWd <- getwd()
setwd('/mrc-bsu/scratch/zmx21/UKB_Data/blood_pressure_test/')
genotype_data <- rbgen::bgen.load(filename = 'sbp_chr10.bgen',rsids = manualResults$rsid)
setwd(curWd)

# perfectImputation <- (genotype_data$data[,,1] == 1 | genotype_data$data[,,1] == 0)
# ratioPerfectImputed <- apply(perfectImputation,1,function(x) sum(x)/length(x))
# ratioPerfectImputed <- data.frame(rsid=rownames(genotype_data$data),ratio_perfect_imputed=ratioPerfectImputed)
deviation <- pmin(1 - genotype_data$data[,,1],genotype_data$data[,,1]) + pmin(1 - genotype_data$data[,,2],genotype_data$data[,,2]) + pmin(1 - genotype_data$data[,,3],genotype_data$data[,,3])
sumOfDeviation <- apply(deviation,1,function(x) sum(x)/3/length(x))
sumOfDeviation <- data.frame(rsid=rownames(genotype_data$data),sum_of_deviation=sumOfDeviation)

allResults <- dplyr::inner_join(plinkResults,manualResults,by=c('rsid'='rsid')) %>% dplyr::inner_join(sumOfDeviation,by=c('rsid'='rsid'))


