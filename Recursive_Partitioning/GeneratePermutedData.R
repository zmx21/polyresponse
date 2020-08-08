#Generate permuted data indices
path <- '/home/zmx21/bsu_scratch/LDL_Project_Data_Aug2019/Random_Forest/rs12916_rs17238484_rs5909_rs2303152_rs10066707_rs2006760_LDLdirect'
testing_set <- readRDS('~/bsu_scratch/LDL_Project_Data_Aug2019/Genotype_Data/test_set.rds')
training_set <- readRDS('~/bsu_scratch/LDL_Project_Data_Aug2019/Genotype_Data/training_set.rds')

n_samples <- length(testing_set)
n_perm = 1000
rand_perm <- lapply(1:n_perm, function(x) sample(1:n_samples,size = n_samples,replace = F))
saveRDS(rand_perm,file = paste0(path,'/randPerm.rds'))