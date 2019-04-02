####################################################################################
#Creates random forest, dependent on the implementation of decision tree in InteractionTree.R
#Input: As listed in message below (line 70)
#Output: .rds file for each tree.
####################################################################################

source('~/MRC_BSU_Internship_LDL/Recursive_Partitioning/InteractionTree.R')
source('~/MRC_BSU_Internship_LDL/Recursive_Partitioning/LoadDosage.R')
source('~/MRC_BSU_Internship_LDL/Recursive_Partitioning/ExtractSubsample.R')
source('~/MRC_BSU_Internship_LDL/Load_Bgen/LoadBgen.R')
library(pbmcapply)
library(partykit)
library(parallel)
library(RSQLite)
library(dplyr)

CreateRandomForestPermInteraction <- function(data,testing_set,
                                              n_bootstrap,n_features,tree_min_size,outpath,
                                              n_cores,data_path,p_val_thresh){
  
  #Construct rf with same parameters as non-permuted trees
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
      curTree <- readRDS(paste0(data_path,'0.75_',tree_min_size,'_',p_val_thresh,'/tree',i,'.rds'))
      bootstrapIndex <- curTree$bootstrapIndex
      outofbagIndex <- curTree$outofbagIndex
      
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
    #Get Total interaction of random tree
    CalculateTotalInteraction(as.data.frame(testing_set$dosageMatrix),bootstrapPartyTree,testing_set$dosageTarget,
                              testing_set$phenotypes,testing_set$covariates)
  },mc.cores = n_cores,ignore.interactive = T)
  #})
  saveRDS(res,paste0(outpath,'interactions.rds'))
}
#Parse arguments
##First read in the arguments listed at the command line
args=(commandArgs(TRUE))
if(length(args)==0){
  interaction_path <- '~/bsu_scratch/LDL_Project_Data/Interaction_Data/HMGCR_LDL_known.txt'
  p_val_thresh <-  '5e-6'
  targetRS <-  'rs12916,rs17238484,rs5909,rs2303152,rs10066707,rs2006760'
  phenotype <- 'LDLdirect'
  outpath <- '~/bsu_scratch/LDL_Project_Data/Random_Forest/Variable_Perm_Interaction/'
  n_bootstrap <- '1:1'
  n_features <- '0.75'
  min_node_size <- 30000
  n_cores <- 1
  perm_n <- 1
  data_path <- '~/bsu_scratch/LDL_Project_Data/Random_Forest/rs12916_rs17238484_rs5909_rs2303152_rs10066707_rs2006760_LDLdirect/'
  targetRS <- unlist(strsplit(targetRS,split = ','))
  outpath <- paste0(outpath,paste(targetRS,collapse = '_'),'_',phenotype,'/')
  system(paste0('mkdir -p ',outpath))
  system(paste0('mkdir -p ',outpath,'/dosage_matrix/'))
  suffix <- paste(n_features,min_node_size,p_val_thresh,sep = '_')
  system(paste0('mkdir -p ',outpath,suffix,'/'))
  system(paste0('mkdir -p ',outpath,suffix,'/perm',perm_n))
  
  
}else if (length(args) != 11){
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
       '#11: Data Path',
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
           "#10: Perm Number",
           "#11: Data Path")
  variables <- c('interaction_path','p_val_thresh','targetRS','phenotype','outpath','n_bootstrap','n_features','min_node_size','n_cores','perm_n','data_path')
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
  data <- readRDS(paste0(data_path,'data_p_',as.numeric(p_val_thresh),'.rds'))
  data$dosageMatrix <- readRDS(paste0(outpath,'/dosage_matrix/','dosage_matrix_perm_',as.numeric(perm_n),'_p_',as.numeric(p_val_thresh),'.rds'))
}else{
  print('creating dosage matrix with random SNPs')
  data <- readRDS(paste0(data_path,'data_p_',as.numeric(p_val_thresh),'.rds'))
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
  path <-  '~/bsu_scratch/LDL_Project_Data/Genotype_Data/'
  bgen_file_prefix <- 'ukb_imp_chr#_HRConly'
  
  randDosageMatrix <- t(as.matrix(LoadBgen(path,bgen_file_prefix,rand_pred)))[rownames(data$dosageMatrix),]
  saveRDS(randDosageMatrix,paste0(outpath,'/dosage_matrix/','dosage_matrix_perm_',as.numeric(perm_n),'_p_',as.numeric(p_val_thresh),'.rds'))
  data$dosageMatrix <- randDosageMatrix
}
#Extract training set.
training_testing_set <- ExtractSubSample(data,
                                         readRDS('~/bsu_scratch/LDL_Project_Data/Genotype_Data/training_set.rds'),
                                         readRDS('~/bsu_scratch/LDL_Project_Data/Genotype_Data/test_set.rds'))
training_set <- training_testing_set$bootstrap
testing_set <- training_testing_set$outofbag
#Create random forest
CreateRandomForestPermInteraction(data = training_set,
                                  testing_set = testing_set,
                                  n_bootstrap = n_bootstrap,
                                  n_features = n_features,
                                  tree_min_size = as.numeric(min_node_size),
                                  outpath = paste0(outpath,suffix,'/perm',perm_n,'/'),
                                  n_cores = as.numeric(n_cores),
                                  data_path = data_path,
                                  p_val_thresh = p_val_thresh)
