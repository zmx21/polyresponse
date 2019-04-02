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

GetTotalInteraction <- function(testing_set,n_bootstrap,outpath,n_cores){
  system(paste0('mkdir -p ',outpath,'tree_interactions/'))
  #Get all trees
  if(grepl(':',n_bootstrap)){
    n_bootstrap <- unlist(strsplit(n_bootstrap,':'))
    chunks <- as.numeric(n_bootstrap[1]):as.numeric(n_bootstrap[2])
  }else{
    chunks <- 1:as.numeric(n_bootstrap)
  }
  res <- pbmclapply(chunks,function(i) {
  #res <- lapply(chunks,function(i){
    curTree <- readRDS(paste0(outpath,'tree',i,'.rds'))
    #Get Total interaction of current tree
    CalculateTotalInteraction(as.data.frame(testing_set$dosageMatrix),curTree$bootstrapPartyTree,testing_set$dosageTarget,
                              testing_set$phenotypes,testing_set$covariates)
  },mc.cores = n_cores,ignore.interactive = T)
  #})
  saveRDS(res,paste0(outpath,'tree_interactions/','interactions.rds'))
}

#Parse arguments
##First read in the arguments listed at the command line
args=(commandArgs(TRUE))
if(length(args)==0){
  interaction_path <- '~/bsu_scratch/LDL_Project_Data/Interaction_Data/HMGCR_LDL_known.txt'
  p_val_thresh <-  '5e-6'
  targetRS <-  'rs12916,rs17238484,rs5909,rs2303152,rs10066707,rs2006760'
  phenotype <- 'LDLdirect'
  outpath <- '~/bsu_scratch/LDL_Project_Data/Random_Forest/'
  n_bootstrap <- '1:1'
  n_features <- '0.75'
  min_node_size <- 30000
  n_cores <- 16
  targetRS <- unlist(strsplit(targetRS,split = ','))
  outpath <- paste0(outpath,paste(targetRS,collapse = '_'),'_',phenotype,'/')
  suffix <- paste(n_features,min_node_size,p_val_thresh,sep = '_')

  
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
  suffix <- paste(n_features,min_node_size,p_val_thresh,sep = '_')
 
  
}
#Check if data has been loaded already, and already save
data <- readRDS(paste0(outpath,'data_p_',as.numeric(p_val_thresh),'.rds'))

#Extract training set.
training_testing_set <- ExtractSubSample(data,readRDS('~/bsu_scratch/LDL_Project_Data/Genotype_Data/training_set.rds'),
                                         readRDS('~/bsu_scratch/LDL_Project_Data/Genotype_Data/test_set.rds'))
testing_set <- training_testing_set$outofbag
#Create random forest
GetTotalInteraction(testing_set,n_bootstrap,paste0(outpath,suffix,'/'),as.numeric(n_cores))