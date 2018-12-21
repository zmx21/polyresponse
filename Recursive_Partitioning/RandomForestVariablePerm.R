####################################################################################
#Creates random forest, dependent on the implementation of decision tree in InteractionTree.R
#Input: As listed in message below (line 70)
#Output: .rds file for each tree.
####################################################################################

source('~/MRC_BSU_Internship/Recursive_Partitioning/InteractionTree.R')
source('~/MRC_BSU_Internship/Recursive_Partitioning/LoadDosage.R')
source('~/MRC_BSU_Internship/Recursive_Partitioning/ExtractSubsample.R')
source('~/MRC_BSU_Internship/Load_Bgen/LoadBgen.R')
library(pbmcapply)
library(partykit)
library(parallel)
library(RSQLite)
library(dplyr)

#Generate training and testing set prediction, given a genotype matrix
GeneratePrediction <- function(tree,genotypeMatrix,testingSetSamples){
  nodeAssignment <- predict(tree,genotypeMatrix)
  nodeAssignment <- as.vector(nodeAssignment)
  #Calculate treatment effect estimates based on training set
  trainingMainEffect <- unlist(nodeapply(tree,ids = unique(nodeAssignment),FUN = function(n) as.numeric(unlist(strsplit(x = unlist(n$info),split = "  "))[1])))
  
  #Calculate treatment effect estimates based on testing set
  uniqueLeafNodes <- unique(nodeAssignment)
  testingMainEffects <- rep(NA,length(uniqueLeafNodes))
  names(testingMainEffects) <- uniqueLeafNodes
  for(i in 1:length(uniqueLeafNodes)){
    currentIndicator <- nodeAssignment == uniqueLeafNodes[i]
    currentSubset <- lapply(testingSetSamples,function(x) SubsetData(x,currentIndicator))
    testingMainEffects[i] <- FitMainEffectModel(currentSubset$dosageTarget,currentSubset$phenotypes,currentSubset$covariates)$coeff['treatment']
  }
  results <- list(nodeAssignment = nodeAssignment,testingMainEffects = testingMainEffects,trainingMainEffect=trainingMainEffect)
  return(results)
}

#Runs predictions of testing set samples
PredictFromVarPermRF <- function(testingSetSamples,randomForestPath,tree,tree_index){
  #Extract genotype matrix from testing set
  genotypeMatrix <- as.data.frame(testingSetSamples$dosageMatrix)
  #Calculate beta coeff for each sample for the tree
  predResult <- GeneratePrediction(tree,genotypeMatrix,testingSetSamples)
  saveRDS(predResult,paste0(randomForestPath,'/prediction_betas/','tree',tree_index,'.rds'))
}

CreateRandomForestPerm <- function(data,testing_set,sample_size,n_bootstrap,n_features,tree_min_size,outpath,n_cores,perm_n){
  system(paste0('mkdir -p ',outpath,'/prediction_betas/'))
  #vector of all sample ids
  samples <- 1:nrow(data$dosageMatrix)
  
  #Construct random forest
  if(grepl(':',n_bootstrap)){
    n_bootstrap <- unlist(strsplit(n_bootstrap,':'))
    chunks <- as.numeric(n_bootstrap[1]):as.numeric(n_bootstrap[2])
  }else{
    chunks <- 1:as.numeric(n_bootstrap)
  }
  res <- pbmclapply(chunks,function(i) {
    #res <- lapply(chunks,function(i){
    flag <- T
    while(flag){
      bootstrapIndex <- samples[sample(1:length(samples),size = sample_size,replace = T)]
      outofbagIndex <- setdiff(samples,bootstrapIndex)
      
      #Extract the boostrap sample and out of bag sample.
      currentSubsample <- ExtractSubSample(data,bootstrapIndex,outofbagIndex)
      currentBootstrap <- currentSubsample$bootstrap
      currentOutOfBag <- currentSubsample$outofbag
      #Check how many features to choose randomly (if expressed as number,fraction, or default)
      if(suppressWarnings(is.na(as.numeric(n_features)))){
        n_features <- floor(eval(parse(text = paste0(n_features,'(',dim(currentBootstrap$dosageMatrix)[2],')'))))
      }else if(as.numeric(n_features) <= 1){
        n_features <- floor(as.numeric(n_features) * dim(currentBootstrap$dosageMatrix)[2])
      }else{
        n_features <- as.numeric(n_features)
      }
      #Construct current decision tree
      bootstrapTree <- tryCatch({
        ConstructTree(currentBootstrap,tree_min_size,n_features)
      }, error = function(e) {
      })
      bootstrapPartyTree <- tryCatch({
        party(bootstrapTree,as.data.frame(currentBootstrap$dosageMatrix))
      }, error = function(e) {
      })
      if(is.list(bootstrapPartyTree)){
        flag <- F
      }
    }
    #Save prediction results
    PredictFromVarPermRF(testing_set,outpath,bootstrapPartyTree,i)
    
  },mc.cores = n_cores,ignore.interactive = T)
  #})
}
#Parse arguments
##First read in the arguments listed at the command line
args=(commandArgs(TRUE))
if(length(args)==0){
  interaction_path <- '~/parsed_interaction/CACNA1D_sbp.txt'
  p_val_thresh <-  '1e-5'
  targetRS <-  'rs3821843,rs7340705,rs113210396,rs312487,rs11719824,rs3774530,rs3821856'
  phenotype <- 'sbp'
  outpath <- '~/bsu_scratch/Random_Forest/Variable_Perm/'
  n_bootstrap <- '1:1'
  n_features <- '0.75'
  min_node_size <- 30000
  n_cores <- 1
  perm_n <- 1
  targetRS <- unlist(strsplit(targetRS,split = ','))
  outpath <- paste0(outpath,paste(targetRS,collapse = '_'),'_',phenotype,'/')
  system(paste0('mkdir -p ',outpath))
  system(paste0('mkdir -p ',outpath,'/dosage_matrix/'))
  suffix <- paste(n_features,min_node_size,p_val_thresh,sep = '_')
  system(paste0('mkdir -p ',outpath,suffix,'/'))
  system(paste0('mkdir -p ',outpath,suffix,'/perm',perm_n))
  
  
}else if (length(args) != 10){
  stop("You need to supply:\n",
       "# 1: Interaction result path\n",
       '# 2: Interaction P-value threshold\n',
       "# 3: Target rsid\n",
       "# 4: Phenotype\n",
       "# 5: RF output path\n",
       "# 6: Number of bootstraps\n",
       "# 7: Number of random features\n",
       "# 8: Minimum node size\n",
       "# 9: Number of cores",
       "#10: Perm Number",
       "Exiting...", call.=FALSE)
}else{
  print('All Arguments Supplied')
  str <- c("# 1: Interaction result path:",
           '# 2: Interaction P-value threshold:',
           "# 3: Target rsid:",
           "# 4: Phenotype:",
           "# 5: RF output path:",
           "# 6: Number of bootstraps:",
           "# 7: Number of random features:",
           "# 8: Minimum node size:",
           "# 9: Number of cores",
           "#10: Perm Number")
  variables <- c('interaction_path','p_val_thresh','targetRS','phenotype','outpath','n_bootstrap','n_features','min_node_size','n_cores','perm_n')
  for(i in 1:length(args)){
    eval(parse(text=paste0(variables[i],'=',"'",args[[i]],"'")))
    print(paste0(str[i],"   ",args[i]))
  }
  #Split rsid based on comma
  targetRS <- unlist(strsplit(targetRS,split = ','))
  
  outpath <- paste0(outpath,paste(targetRS,collapse = '_'),'_',phenotype,'/')
  system(paste0('mkdir -p ',outpath))
  system(paste0('mkdir -p ',outpath,'/dosage_matrix/'))
  suffix <- paste(n_features,min_node_size,p_val_thresh,sep = '_')
  system(paste0('mkdir -p ',outpath,suffix,'/'))
  system(paste0('mkdir -p ',outpath,suffix,'/perm',perm_n))
}
#Check if data has been loaded already, and already saved
if(paste0('dosage_matrix_perm_',as.numeric(perm_n),'_p_',as.numeric(p_val_thresh),'.rds') %in% dir(paste0(outpath,'/dosage_matrix/'))){
  print('loading existing dosage matrix')
  data <- readRDS(paste0('/home/zmx21/bsu_scratch/Random_Forest/rs3821843_rs7340705_rs113210396_rs312487_rs11719824_rs3774530_rs3821856_sbp/','data_p_',as.numeric(p_val_thresh),'.rds'))
  data$dosageMatrix <- readRDS(paste0(outpath,'/dosage_matrix/','dosage_matrix_perm_',as.numeric(perm_n),'_p_',as.numeric(p_val_thresh),'.rds'))
}else{
  print('creating dosage matrix with random SNPs')
  data <- readRDS(paste0('/home/zmx21/bsu_scratch/Random_Forest/rs3821843_rs7340705_rs113210396_rs312487_rs11719824_rs3774530_rs3821856_sbp/','data_p_',as.numeric(p_val_thresh),'.rds'))
  true_predictors <- colnames(data$dosageMatrix)
  n_predictors <- length(true_predictors)
  interaction_results <- data.table::fread(interaction_path)
  #choose 2000 random SNPs
  interaction_results <- interaction_results[sample(1:nrow(interaction_results),size = 2000),]
  #Connect to rsid annotation database. 
  anno_sql_name<- "all_snp_stats.sqlite"
  path <- '~/bsu_scratch/SQL/'
  curwd <- getwd()
  setwd(path)
  anno_con <- RSQLite::dbConnect(SQLite(), dbname = anno_sql_name)
  anno_db <- dplyr::tbl(anno_con,'all_snp_stats')
  #Get all rsids meeting criteria
  true_anno <- anno_db %>% dplyr::filter(rsid %in% true_predictors) %>% collect() %>% dplyr::select(rsid,minor_allele_frequency)
  anno <- anno_db %>% dplyr::filter(rsid %in% interaction_results$rsid) %>% collect()
  dbDisconnect(anno_con)
  setwd(curwd)
  interaction_results <- interaction_results %>% dplyr::left_join(anno,by=c('rsid'='rsid')) %>% dplyr::filter(as.numeric(info) > 0.5)
  
  #Find matching pred with similar MAF
  rand_pred <- c()
  for(i in 1:nrow(true_anno)){
    currentRS <- true_anno$rsid[i]
    currentMAF <- as.numeric(true_anno$minor_allele_frequency[i])
    closestRS <- interaction_results$rsid[which.min(abs(currentMAF - as.numeric(interaction_results$minor_allele_frequency)))]
    rand_pred <- c(rand_pred,closestRS)
  }
  path <-  '/mrc-bsu/scratch/zmx21/UKB_Data/'
  bgen_file_prefix <- 'ukb_imp_chr#_HRConly'
  
  randDosageMatrix <- t(as.matrix(LoadBgen(path,bgen_file_prefix,rand_pred)))[rownames(data$dosageMatrix),]
  saveRDS(randDosageMatrix,paste0(outpath,'/dosage_matrix/','dosage_matrix_perm_',as.numeric(perm_n),'_p_',as.numeric(p_val_thresh),'.rds'))
  data$dosageMatrix <- randDosageMatrix
}
#Extract training set.
training_testing_set <- ExtractSubSample(data,readRDS('~/bsu_scratch/UKB_Data/training_set.rds'),readRDS('~/bsu_scratch/UKB_Data/test_set.rds'))
training_set <- training_testing_set$bootstrap
testing_set <- training_testing_set$outofbag
#Create random forest
CreateRandomForestPerm(training_set,testing_set,floor(nrow(training_set$dosageMatrix)*(2/3)),n_bootstrap,n_features,as.numeric(min_node_size),paste0(outpath,suffix,'/perm',perm_n,'/'),as.numeric(n_cores))
