CompareInteractions <- function(targetRS,suffix,n_perm){
  print(suffix)
  #Get sum of interactions from true tree
  trueTreeInteractions <- readRDS(paste0('~/bsu_scratch/LDL_Project_Data/Random_Forest/',targetRS,suffix,'/tree_interactions/interactions.rds'))
  #trueMeanSumOfInteractions <- mean(sapply(trueTreeInteractions,function(x) sum(x)))
  
  #Get sum of interactions from tree with random predictors
  permTreeInteractions <- lapply(1:2000,function(i) readRDS(paste0('~/bsu_scratch/LDL_Project_Data/Random_Forest/',
                                                             targetRS,suffix,'interactions_pheno_perm/interactions_tree',i,'.rds')))
  #permMeanSumOfInteractions <- sapply(permTreeInteractions,function(x) mean(sapply(x,function(y) sum(y))))
  return(list(true=trueTreeInteractions,perm=permTreeInteractions))
}


targetRS <- 'rs12916_rs17238484_rs5909_rs2303152_rs10066707_rs2006760_LDLdirect/'
node_size <- seq(10000,40000,10000)
thresh <- c('5e-6','1e-5','3e-5','5e-5')
#thresh <- c('5e-6','3e-5')
comb <- expand.grid(node_size,thresh)
colnames(comb) <- c('node_size','thresh')
n_perm <- 50

# res <- lapply(1:nrow(comb),function(i) CompareInteractions(targetRS,paste0('0.75_',comb$node_size[i],'_',as.character(comb$thresh[i]),'/'),n_perm))
# saveRDS(res,file='~/bsu_scratch/LDL_Project_Data/Random_Forest/rs12916_rs17238484_rs5909_rs2303152_rs10066707_rs2006760_LDLdirect/interaction_pheno_perm.rds')

# res <- readRDS('~/bsu_scratch/LDL_Project_Data/Random_Forest/rs12916_rs17238484_rs5909_rs2303152_rs10066707_rs2006760_LDLdirect/interaction_pheno_perm.rds')
n_comb <- 16
n_tree=2000


rel_pop_diff <- data.frame(thresh=numeric(),node_size=numeric(),mean_diff=numeric(),sd_diff=numeric(),low_CI=numeric(),high_CI=numeric())
abs_diff <- data.frame(thresh=numeric(),node_size=numeric(),mean_diff=numeric(),sd_diff=numeric(),low_CI=numeric(),high_CI=numeric())
perm_true_diff<- data.frame(thresh=numeric(),node_size=numeric(),mean_diff=numeric(),sd_diff=numeric(),low_CI=numeric(),high_CI=numeric())
for(i in 1:n_comb){
  curCombPerm <- res[[i]]$perm
  curCombTrue <- res[[i]]$true
  
  trueTreeTrainTestDiff <- 0
  truePermDiff <- rep(0,n_perm)
  permTrueDiff <- rep(0,n_perm)
  for(j in 1:n_tree){
    curTreePerm <- curCombPerm[[j]]
    curTreeTrue <- curCombTrue[[j]]
    
    
    curTrueTreeTrainTestDiff <- sum(curTreeTrue$train - curTreeTrue$test)
    curPermTreeTrainTestDiff <- sapply(curTreePerm,function(x) sum(curTreeTrue$train - x))
   
    # curSumTreePerm <- sapply(curTreePerm,function(x) sum(x))
    # curTrueTreeTrainTestDiff <- sum(curTreeTrue$train) - sum(curTreeTrue$test)
    # curPermTreeTrainTestDiff <- sum(curTreeTrue$train) - sumTreePerm
    
    curTruePermDiff <- curTrueTreeTrainTestDiff - curPermTreeTrainTestDiff
    
    trueTreeTrainTestDiff <- trueTreeTrainTestDiff + curTrueTreeTrainTestDiff
    truePermDiff <- truePermDiff + curTruePermDiff
    
    
    curPermTrueDiff <- sapply(curTreePerm,function(x) sum(curTreeTrue$test) - sum(x))
    permTrueDiff <- permTrueDiff + curPermTrueDiff
  }
  
  popRelDiff <- truePermDiff / trueTreeTrainTestDiff
  absDiff <- truePermDiff / n_tree
  rel_pop_diff <- rbind(rel_pop_diff,data.frame(thresh=comb$thresh[i],
                                                node_size=comb$node_size[i],
                                                mean_diff=mean(popRelDiff),
                                                sd_diff=sd(popRelDiff),
                                                low_CI=quantile(popRelDiff,probs = c(0.025)),
                                                high_CI=quantile(popRelDiff,probs = c(0.975))))
  
  abs_diff <- rbind(abs_diff,data.frame(thresh=comb$thresh[i],
                                                node_size=comb$node_size[i],
                                                mean_diff=mean(absDiff),
                                                sd_diff=sd(absDiff),
                                                low_CI=quantile(absDiff,probs = c(0.025)),
                                                high_CI=quantile(absDiff,probs = c(0.975))))
  
  perm_true_diff <- rbind(perm_true_diff,data.frame(thresh=comb$thresh[i],
                                              node_size=comb$node_size[i],
                                              mean_diff=mean(permTrueDiff),
                                              sd_diff=sd(permTrueDiff),
                                              low_CI=quantile(permTrueDiff,probs = c(0.025)),
                                              high_CI=quantile(permTrueDiff,probs = c(0.975))))
}


test_train_int_diff_nonperm <- lapply(res,function(x) sapply(x$true,function(y) sum(y$train) - sum(y$test)))
test_int_perm <- lapply(res,function(x) lapply(x$perm,function(y) sapply(y,function(z) sum(z))))

test_train_int_diff_perm_nonperm <- mapply(function(x,y) mapply(function(a,b) mapply(function(b,c) b - c,a,b,SIMPLIFY = F),x,y,SIMPLIFY=F),test_train_int_diff_nonperm,test_int_perm,SIMPLIFY = F)

