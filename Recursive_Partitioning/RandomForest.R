source('~/MRC_BSU_Internship/Recursive_Partitioning/InteractionTree.R')
source('~/MRC_BSU_Internship/Recursive_Partitioning/LoadDosage.R')
source('~/MRC_BSU_Internship/Recursive_Partitioning/ExtractSubsample.R')
library(pbmcapply)
library(partykit)
library(parallel)

CreateRandomForest <- function(data,sample_size,n_bootstrap,n_features,tree_min_size,outpath,n_cores){
  #vector of all sample ids
  samples <- 1:nrow(data$dosageMatrix)

  #Construct random forest
  if(grepl(':',n_bootstrap)){
    n_bootstrap <- unlist(strsplit(n_bootstrap,':'))
    chunks <- as.numeric(n_bootstrap[1]):as.numeric(n_bootstrap[2])
  }else{
    chunks <- 1:as.numeric(n_bootstrap)
  }
  trash <- pbmclapply(chunks,function(i) {
    #create list of boostrap samples to construct tree on
    bootstrapIndex <- samples[sample(1:length(samples),size = sample_size,replace = T)]
    outofbagIndex <- setdiff(samples,bootstrapIndex)
    
    #Extract the boostrap sample and out of bag sample.
    currentSubsample <- ExtractSubSample(data,bootstrapIndex,outofbagIndex)
    currentBootstrap <- currentSubsample$bootstrap
    currentOutOfBag <- currentSubsample$outofbag
    if(suppressWarnings(is.na(as.numeric(n_features)))){
      n_features <- floor(eval(parse(text = paste0(n_features,'(',dim(currentBootstrap$dosageMatrix)[2],')'))))
    }else if(as.numeric(n_features) <= 1){
      n_features <- floor(as.numeric(n_features) * dim(currentBootstrap$dosageMatrix)[2])
    }else{
      n_features <- as.numeric(n_features)
    }
    #Construct decision tree 
    bootstrapTree <- ConstructTree(currentBootstrap,tree_min_size,n_features)
    bootstrapPartyTree <- party(bootstrapTree,as.data.frame(currentBootstrap$dosageMatrix))
    bootstrapTreeObj <- list(bootstrapPartyTree = bootstrapPartyTree,bootstrapIndex=bootstrapIndex,outofbagIndex=outofbagIndex)
    saveRDS(bootstrapTreeObj,paste0(outpath,'tree',i,'.rds'))
  },mc.cores = n_cores,ignore.interactive = T)
}
#Parse arguments
##First read in the arguments listed at the command line
args=(commandArgs(TRUE))
if(length(args)==0){
  interaction_path <- '~/parsed_interaction/CACNA1D_sbp.txt'
  p_val_thresh <-  '1e-5'
  targetRS <-  'rs3821843,rs7340705,rs113210396,rs312487,rs11719824,rs3774530,rs3821856'
  phenotype <- 'sbp'
  outpath <- '~/bsu_scratch/Random_Forest/'
  n_bootstrap <- '1:500'
  n_features <- '0.75'
  min_node_size <- 40000
  n_cores <- 16
  targetRS <- unlist(strsplit(targetRS,split = ','))
  outpath <- paste0(outpath,paste(targetRS,collapse = '_'),'_',phenotype,'/')
  system(paste0('mkdir -p ',outpath))
  suffix <- paste(n_features,min_node_size,p_val_thresh,sep = '_')
  system(paste0('mkdir -p ',outpath,suffix,'/'))
  
}else if (length(args) != 9){
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
           "# 9: Number of cores")
  variables <- c('interaction_path','p_val_thresh','targetRS','phenotype','outpath','n_bootstrap','n_features','min_node_size','n_cores')
  for(i in 1:length(args)){
    eval(parse(text=paste0(variables[i],'=',"'",args[[i]],"'")))
    print(paste0(str[i],"   ",args[i]))
  }
  #Split rsid based on comma
  targetRS <- unlist(strsplit(targetRS,split = ','))
  
  outpath <- paste0(outpath,paste(targetRS,collapse = '_'),'_',phenotype,'/')
  system(paste0('mkdir -p ',outpath))
  suffix <- paste(n_features,min_node_size,p_val_thresh,sep = '_')
  system(paste0('mkdir -p ',outpath,suffix,'/'))
  
}
if(paste0('data_p_',as.numeric(p_val_thresh),'.rds') %in% dir(outpath)){
  data <- readRDS(paste0(outpath,'data_p_',as.numeric(p_val_thresh),'.rds'))
}else{
  data <- LoadDosage(as.numeric(p_val_thresh),interaction_path,phenotype,0.3,0.05,0.5,targetRS)
  saveRDS(data,file = paste0(outpath,'data_p_',as.numeric(p_val_thresh),'.rds'))
}
training_testing_set <- ExtractSubSample(data,readRDS('~/bsu_scratch/UKB_Data/training_set.rds'),readRDS('~/bsu_scratch/UKB_Data/test_set.rds'))
training_set <- training_testing_set$bootstrap
CreateRandomForest(training_set,floor(nrow(training_set$dosageMatrix)*(2/3)),n_bootstrap,n_features,as.numeric(min_node_size),paste0(outpath,suffix,'/'),as.numeric(n_cores))
