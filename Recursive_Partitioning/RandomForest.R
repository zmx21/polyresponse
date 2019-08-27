####################################################################################
#Creates random forest, dependent on the implementation of decision tree in InteractionTree.R
#Input: As listed in message below (line 70)
#Output: .rds file for each tree.
####################################################################################

source('~/MRC_BSU_Internship_LDL/Recursive_Partitioning/InteractionTree.R')
source('~/MRC_BSU_Internship_LDL/Recursive_Partitioning/LoadDosage.R')
source('~/MRC_BSU_Internship_LDL/Recursive_Partitioning/ExtractSubsample.R')
library(pbmcapply)
library(partykit)
library(parallel)

CreateRandomForest <- function(data,sample_size,n_bootstrap,n_features,tree_min_size,outpath,bootStrapIndex,n_cores){
  #Construct random forest
  if(grepl(':',n_bootstrap)){
    n_bootstrap <- unlist(strsplit(n_bootstrap,':'))
    chunks <- as.numeric(n_bootstrap[1]):as.numeric(n_bootstrap[2])
  }else{
    chunks <- 1:as.numeric(n_bootstrap)
  }
  res <- pbmclapply(chunks,function(i) {
  #res <- lapply(chunks,function(i){
    #extract list of boostrap samples (sample with replacement) to construct current tree on
    curBootstrapIndex <- bootStrapIndex[[i]]
    curOutOfbagIndex <- setdiff(1:nrow(data$dosageMatrix),curBootstrapIndex)
    flag <- T
    while(flag){
      #Extract the boostrap sample and out of bag sample.
      currentSubsample <- ExtractSubSample(data,curBootstrapIndex,curOutOfbagIndex)
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
    bootstrapTreeObj <- list(bootstrapPartyTree = bootstrapPartyTree,bootstrapIndex=curBootstrapIndex,outofbagIndex=curOutOfbagIndex)
    #Save tree, and bootstrap indices
    saveRDS(bootstrapTreeObj,paste0(outpath,'tree',i,'.rds'))
  },mc.cores = n_cores,ignore.interactive = T)
  #})
}
#Parse arguments
##First read in the arguments listed at the command line
args=(commandArgs(TRUE))
if(length(args)==0){
  interaction_path <- '~/bsu_scratch/LDL_Project_Data_Aug2019/Interaction_Data/HMGCR_LDL_known.txt'
  p_val_thresh <-  '5e-6'
  targetRS <-  'rs12916,rs17238484,rs5909,rs2303152,rs10066707,rs2006760'
  phenotype <- 'LDLdirect'
  outpath <- '~/bsu_scratch/LDL_Project_Data_Aug2019/Random_Forest/'
  n_bootstrap <- '1:2000'
  n_features <- '0.75'
  min_node_size <- 40000
  n_cores <- 16
  MAF <- '5e-2'
  targetRS <- unlist(strsplit(targetRS,split = ','))
  outpath <- paste0(outpath,paste(targetRS,collapse = '_'),'_',phenotype,'/')
  system(paste0('mkdir -p ',outpath))
  suffix <- paste(n_features,min_node_size,p_val_thresh,MAF,sep = '_')
  system(paste0('mkdir -p ',outpath,suffix,'/'))
  
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
       "# 10: MAF",
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
           '#10: MAF')
  variables <- c('interaction_path','p_val_thresh','targetRS','phenotype','outpath','n_bootstrap','n_features','min_node_size','n_cores','MAF')
  for(i in 1:length(args)){
    eval(parse(text=paste0(variables[i],'=',"'",args[[i]],"'")))
    print(paste0(str[i],"   ",args[i]))
  }
  #Split rsid based on comma
  targetRS <- unlist(strsplit(targetRS,split = ','))
  
  outpath <- paste0(outpath,paste(targetRS,collapse = '_'),'_',phenotype,'/')
  system(paste0('mkdir -p ',outpath))
  suffix <- paste(n_features,min_node_size,p_val_thresh,MAF,sep = '_')
  system(paste0('mkdir -p ',outpath,suffix,'/'))
  
}
CheckData <- function(outpath,p_val_thresh,interaction_path,phenotype,targetRS,MAF,trainingSize,percBootstrap,nBootstraps){
  #Check if data has been loaded already, and already saved
  if(paste0('data_p_',as.numeric(p_val_thresh),'_maf_',MAF,'.rds') %in% dir(outpath)){
    data <- readRDS(paste0(outpath,'data_p_',as.numeric(p_val_thresh),'_maf_',MAF,'.rds'))
    bootstraps <- readRDS(paste0(outpath,'bootstrap_indices.rds'))
  }else{
    if(grepl(':',nBootstraps)){
      nBootstraps <- unlist(strsplit(nBootstraps,':'))
      nBootstraps <- as.numeric(nBootstraps[2])
    }else{
      nBootstraps <- as.numeric(nBootstraps)
    }
    
    #Load Pheno, target geno, variant geno, and covaraiates
    data <- LoadDosage(as.numeric(p_val_thresh),interaction_path,phenotype,0.3,as.numeric(MAF),0.5,targetRS)
    #Create bootstrap indices
    if(!'bootstrap_indices.rds' %in% dir(outpath)){
      bootstraps <- lapply(1:nBootstraps,function(x) sample(1:trainingSize,size = floor(trainingSize * percBootstrap),replace = T))
      saveRDS(bootstraps,paste0(outpath,'bootstrap_indices.rds'))
    }else{
      bootstraps <- readRDS(paste0(outpath,'bootstrap_indices.rds'))
    }
    saveRDS(data,file = paste0(outpath,'data_p_',as.numeric(p_val_thresh),'_maf_',MAF,'.rds'))
  }
  return(list(data=data,bootstraps=bootstraps))
}
percBootstrap <- 3/4
#Extract training set and testing set indices
traininingSet <- readRDS('~/bsu_scratch/LDL_Project_Data/Genotype_Data/training_set.rds')
testingSet <- readRDS('~/bsu_scratch/LDL_Project_Data/Genotype_Data/test_set.rds')
#Save data, along with bootstraps
dataAndBootstraps <- CheckData(outpath,p_val_thresh,interaction_path,phenotype,targetRS,MAF,length(traininingSet),percBootstrap,n_bootstrap)
data <- dataAndBootstraps$data
bootstraps <- dataAndBootstraps$bootstraps
#Extract training set data
trainingTestingSetData <- ExtractSubSample(data,traininingSet,testingSet)
trainingSetData <- trainingTestingSetData$bootstrap
#Create random forest
CreateRandomForest(trainingSetData,floor(length(traininingSet)*percBootstrap),n_bootstrap,n_features,as.numeric(min_node_size),paste0(outpath,suffix,'/'),bootstraps,as.numeric(n_cores))
