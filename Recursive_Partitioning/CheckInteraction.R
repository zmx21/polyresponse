set.seed(888)
source('~/MRC_BSU_Internship_LDL/Recursive_Partitioning/ExtractSubsample.R')
#Read in genotype and phenotype, indep set of interaction variants
data <- readRDS('~/bsu_scratch/LDL_Project_Data_Aug2019/Random_Forest/rs12916_rs17238484_rs5909_rs2303152_rs10066707_rs2006760_LDLdirect/data_p_9e-06_maf_5e-2.rds')
testingSet <- readRDS('~/bsu_scratch/LDL_Project_Data_Aug2019/Genotype_Data/test_set.rds') 
trainingSet <- readRDS('~/bsu_scratch/LDL_Project_Data_Aug2019/Genotype_Data/training_set.rds') 

#Obtain top snp and a weak snp
interactionData <- data.table::fread('~/bsu_scratch/LDL_Project_Data_Aug2019/Interaction_Data/HMGCR_LDL_known.txt')
weakSNP <- dplyr::filter(interactionData,p_int > 0.9)
weakSNP <- weakSNP[sample(1:nrow(weakSNP),100),]
topSNP <- dplyr::filter(interactionData,rsid %in% colnames(data$dosageMatrix))
topSNP <- topSNP$rsid[which.min(topSNP$p_int)]


#Training set interaction effect
summary(lm(formula = 'LDL ~ HMGCR*topSNP + sex + ages + bmi + PC1 + PC2 +PC3 + PC4 + PC5',data = cbind(data.frame(LDL = data$phenotypes[trainingSet],HMGCR = data$dosageTarget[trainingSet],topSNP = data$dosageMatrix[trainingSet,topSNP]),data$covariates[trainingSet,])))


#Testing set interactions effect
summary(lm(formula = 'LDL ~ HMGCR*topSNP + sex + ages + bmi + PC1 + PC2 +PC3 + PC4 + PC5',data = cbind(data.frame(LDL = data$phenotypes[testingSet],HMGCR = data$dosageTarget[testingSet],topSNP = data$dosageMatrix[testingSet,topSNP]),data$covariates[testingSet,])))

#Segregate samples based on alelles of top SNP
trainingSetData <- ExtractSubSample(data,trainingSet,testingSet)$outofbag
homoRef <- trainingSetData$dosageMatrix[,topSNP]==0
homoRefFit <- summary(lm(formula = 'LDL ~ HMGCR + sex + ages + bmi + PC1 + PC2 +PC3 + PC4 + PC5',data = cbind(data.frame(LDL = trainingSetData$phenotypes[homoRef],HMGCR = trainingSetData$dosageTarget[homoRef]),trainingSetData$covariates[homoRef,])))
nonHomoRefFit <- summary(lm(formula = 'LDL ~ HMGCR + sex + ages + bmi + PC1 + PC2 +PC3 + PC4 + PC5',data = cbind(data.frame(LDL = trainingSetData$phenotypes[!homoRef],HMGCR = trainingSetData$dosageTarget[!homoRef]),trainingSetData$covariates[!homoRef,])))
trueDiff <- abs(homoRefFit$coefficients['HMGCR','Estimate'] - nonHomoRefFit$coefficients['HMGCR','Estimate'])

#Segregate Sample Randomly (only phenotype)
nPerm <- 100
diffRandPheno <- rep(NaN,nPerm)
permFit1 <- list()
permFit2 <- list()
for(i in 1:nPerm){
  randSeg <- sample(homoRef,size = length(homoRef),replace = F)
  permFit1[[i]] <- summary(lm(formula = 'LDL ~ HMGCR + sex + ages + bmi + PC1 + PC2 +PC3 + PC4 + PC5',data = cbind(data.frame(LDL = trainingSetData$phenotypes[randSeg],HMGCR = trainingSetData$dosageTarget[homoRef]),trainingSetData$covariates[homoRef,])))
  permFit2[[i]] <- summary(lm(formula = 'LDL ~ HMGCR + sex + ages + bmi + PC1 + PC2 +PC3 + PC4 + PC5',data = cbind(data.frame(LDL = trainingSetData$phenotypes[!randSeg],HMGCR = trainingSetData$dosageTarget[!homoRef]),trainingSetData$covariates[!homoRef,])))
  diffRandPheno[i] <- abs(permFit1[[i]]$coefficients['HMGCR','Estimate'] - permFit2[[i]]$coefficients['HMGCR','Estimate'])
}

#Very weak SNP
source('~/MRC_BSU_Internship_LDL/Load_Phenotype/Load_Phenotype.R')
source('~/MRC_BSU_Internship_LDL/Load_Bgen/LoadBgen.R')
#Load phenotype information and covariates
path <-  '~/bsu_scratch/LDL_Project_Data/Genotype_Data/'
bgen_file_prefix <- 'ukb_imp_chr#_HRConly'
sample_file_prefix <- 'ukbb_LDL_metadata_with_PC'
cov_names <- c('sex','ages','bmi')
eur_only <- 1

#Load phenotype and samples to keep
phenotypesAndCov <- LoadPhenotype(path,sample_file_prefix,'LDLdirect',c(cov_names,'PC1','PC2','PC3','PC4','PC5'),
                                  eur_only=1,med=1)
samplesToKeep <- phenotypesAndCov$samplesToKeep

#Segregate samples based on alelles of weak SNPs
weakDiff <- rep(NaN,nrow(weakSNP))
for(i in 1:nrow(weakSNP)){
  weakDosage <- LoadBgen(path,bgen_file_prefix,weakSNP$rsid[i])
  weakDosage <- weakDosage[,samplesToKeep][trainingSet]
  homoRefWeak <- weakDosage==0
  homoRefWeakFit <- summary(lm(formula = 'LDL ~ HMGCR + sex + ages + bmi + PC1 + PC2 +PC3 + PC4 + PC5',data = cbind(data.frame(LDL = trainingSetData$phenotypes[homoRefWeak],HMGCR = trainingSetData$dosageTarget[homoRefWeak]),trainingSetData$covariates[homoRefWeak,])))
  nonHomoRefWeakFit <- summary(lm(formula = 'LDL ~ HMGCR + sex + ages + bmi + PC1 + PC2 +PC3 + PC4 + PC5',data = cbind(data.frame(LDL = trainingSetData$phenotypes[!homoRefWeak],HMGCR = trainingSetData$dosageTarget[!homoRefWeak]),trainingSetData$covariates[!homoRefWeak,])))
  weakDiff[i] <- abs(homoRefWeakFit$coefficients['HMGCR','Estimate'] - nonHomoRefWeakFit$coefficients['HMGCR','Estimate'])
}
