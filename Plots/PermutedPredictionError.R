library(pbmcapply)
library(dbplyr)
#Plot comparison of root squared error between true testing set and permuted testing set
CalcPredErrorDiffRel <- function(path){
  print(path)
  nonPermError <- readRDS(paste0(path,'beta_err.rds'))
  permErrorFiles <- dir(path,pattern = 'beta_err_pheno_perm')
  if(length(permErrorFiles) > 1){
    chunks <- lapply(permErrorFiles,function(x) strsplit(x,'_')[[1]])
    chunks <- lapply(chunks,function(x) sapply(x,function(y) gsub(x = y,pattern = '.rds',replacement = '',fixed = T)))
    chunks <- lapply(chunks,function(x) as.numeric(x)[!is.na(as.numeric(x))])
    chunks <- lapply(chunks,function(x) x[1]:x[2])
  }else{
    chunks <- list(c(1:2000))
    permError <- readRDS(paste0(path,'beta_err_pheno_perm.rds'))
  }
  relDiff <- rep(0,1000)
  pastPermFile <- 0

  numSubGroup <- 0
  for(i in 1:length(nonPermError)){
    curNonPerm <- nonPermError[[i]]
    if(length(permErrorFiles) > 1){
      curPermFile <- which(sapply(chunks,function(x) i%in%x))
      if(curPermFile != pastPermFile){
        permError <- readRDS(paste0(path,permErrorFiles[curPermFile]))
      }
      pastPermFile <- curPermFile
      permIndex <- i - min(chunks[[curPermFile]]) + 1
    }else{
      permIndex <- i
    }
    curPerm <- permError[[permIndex]]
    curDiff <- sapply(curPerm,function(x) sum((curNonPerm - x) / curNonPerm))
    relDiff <- relDiff + curDiff
    numSubGroup <- numSubGroup + length(curNonPerm)
  }
  return(relDiff/numSubGroup)
}
CalcPredErrorDiffRelPop <- function(path){
  print(path)
  nonPermError <- readRDS(paste0(path,'beta_err.rds'))
  permErrorFiles <- dir(path,pattern = 'beta_err_pheno_perm')
  if(length(permErrorFiles) > 1){
    chunks <- lapply(permErrorFiles,function(x) strsplit(x,'_')[[1]])
    chunks <- lapply(chunks,function(x) sapply(x,function(y) gsub(x = y,pattern = '.rds',replacement = '',fixed = T)))
    chunks <- lapply(chunks,function(x) as.numeric(x)[!is.na(as.numeric(x))])
    chunks <- lapply(chunks,function(x) x[1]:x[2])
  }else{
    chunks <- list(c(1:2000))
    permError <- readRDS(paste0(path,'beta_err_pheno_perm.rds'))
  }
  
  relDiff <- rep(0,1000)
  norm <- 0
  
  numSubGroup <- 0
  pastPermFile <- 0
  for(i in 1:length(nonPermError)){
    curNonPerm <- nonPermError[[i]]
    if(length(permErrorFiles) > 1){
      curPermFile <- which(sapply(chunks,function(x) i%in%x))
      if(curPermFile != pastPermFile){
        permError <- readRDS(paste0(path,permErrorFiles[curPermFile]))
      }
      pastPermFile <- curPermFile
      permIndex <- i - min(chunks[[curPermFile]]) + 1
    }else{
      permIndex <- i
    }
    curPerm <- permError[[permIndex]]
    curDiff <- sapply(curPerm,function(x) sum(curNonPerm - x))
    curNorm <- sum(curNonPerm)
    
    relDiff <- relDiff + curDiff
    norm <- norm + curNorm
    numSubGroup <- numSubGroup + length(curNonPerm)
  }
  return(relDiff/norm)
}

CalcPredErrorDiffAbs <- function(path){
  print(path)
  nonPermError <- readRDS(paste0(path,'beta_err.rds'))
  
  permErrorFiles <- dir(path,pattern = 'beta_err_pheno_perm')
  if(length(permErrorFiles) > 1){
    chunks <- lapply(permErrorFiles,function(x) strsplit(x,'_')[[1]])
    chunks <- lapply(chunks,function(x) sapply(x,function(y) gsub(x = y,pattern = '.rds',replacement = '',fixed = T)))
    chunks <- lapply(chunks,function(x) as.numeric(x)[!is.na(as.numeric(x))])
    chunks <- lapply(chunks,function(x) x[1]:x[2])
  }else{
    chunks <- list(c(1:2000))
    permError <- readRDS(paste0(path,'beta_err_pheno_perm.rds'))
  }
  
  diff <- rep(0,1000)
  
  numSubGroup <- 0
  pastPermFile <- 0
  for(i in 1:length(nonPermError)){
    curNonPerm <- nonPermError[[i]]
    if(length(permErrorFiles) > 1){
      curPermFile <- which(sapply(chunks,function(x) i%in%x))
      if(curPermFile != pastPermFile){
        permError <- readRDS(paste0(path,permErrorFiles[curPermFile]))
      }
      pastPermFile <- curPermFile
      permIndex <- i - min(chunks[[curPermFile]]) + 1
    }else{
      permIndex <- i
    }
    curPerm <- permError[[permIndex]]
    curDiff <- sapply(curPerm,function(x) sum((curNonPerm - x)))
    diff <- diff + curDiff
    numSubGroup <- numSubGroup + length(curNonPerm)
  }
  return(diff/numSubGroup)
}

CalcPredErrorAbs <- function(path){
  print(path)
  nonPermError <- readRDS(paste0(path,'beta_err.rds'))
  permErrorFiles <- dir(path,pattern = 'beta_err_pheno_perm')
  if(length(permErrorFiles) > 1){
    chunks <- lapply(permErrorFiles,function(x) strsplit(x,'_')[[1]])
    chunks <- lapply(chunks,function(x) sapply(x,function(y) gsub(x = y,pattern = '.rds',replacement = '',fixed = T)))
    chunks <- lapply(chunks,function(x) as.numeric(x)[!is.na(as.numeric(x))])
    chunks <- lapply(chunks,function(x) x[1]:x[2])
  }else{
    chunks <- list(c(1:2000))
    permError <- readRDS(paste0(path,'beta_err_pheno_perm.rds'))
  }
  nonPermRMSE <- 0
  permRMSE <- rep(0,1000)
  
  numSubGroup <- 0
  pastPermFile <- 0
  for(i in 1:length(nonPermError)){
    curNonPerm <- nonPermError[[i]]
    if(length(permErrorFiles) > 1){
      curPermFile <- which(sapply(chunks,function(x) i%in%x))
      if(curPermFile != pastPermFile){
        permError <- readRDS(paste0(path,permErrorFiles[curPermFile]))
      }
      pastPermFile <- curPermFile
      permIndex <- i - min(chunks[[curPermFile]]) + 1
    }else{
      permIndex <- i
    }
    curPerm <- permError[[permIndex]]
    numSubGroup <- numSubGroup + length(curNonPerm)
    permRMSE <- permRMSE + sapply(curPerm,function(x) sum(x))
    nonPermRMSE <-nonPermRMSE + sum(curNonPerm)
  }
  return(list(permRMSE=permRMSE/numSubGroup,nonPermRMSE=nonPermRMSE/numSubGroup))
}



resultPath <- '~/bsu_scratch/LDL_Project_Data/Random_Forest/rs12916_rs17238484_rs5909_rs2303152_rs10066707_rs2006760_LDLdirect/'
node_size <- c(1000,2500,5000,10000,20000,30000,40000)
thresh <- c('5e-6','1e-5','3e-5','5e-5','7e-5','1e-4')
comb <- expand.grid(node_size,thresh)
colnames(comb) <- c('node_size','thresh')

pred_error_abs <- lapply(1:nrow(comb),function(i) CalcPredErrorAbs(paste0(resultPath,'0.75_',comb$node_size[i],'_',as.character(comb$thresh[i]),'/')))
                                                                   
print('Calc abs diff')
pred_diff_abs <- lapply(1:nrow(comb),function(i) CalcPredErrorDiffAbs(paste0(resultPath,'0.75_',comb$node_size[i],'_',as.character(comb$thresh[i]),'/')))
print('Calc rel diff')
pred_diff_rel <- lapply(1:nrow(comb),function(i) CalcPredErrorDiffRel(paste0(resultPath,'0.75_',comb$node_size[i],'_',as.character(comb$thresh[i]),'/')))
print('Calc rel pop diff')
pred_diff_rel_pop <- lapply(1:nrow(comb),function(i) CalcPredErrorDiffRelPop(paste0(resultPath,'0.75_',comb$node_size[i],'_',as.character(comb$thresh[i]),'/')))


pred_diff_rel_df <- data.frame(p_thresh=numeric(),node_size=numeric(),diff=numeric(),sd=numeric())
for(i in 1:length(pred_diff_rel)){
  curDiff <- as.numeric(pred_diff_rel[[i]])
  pred_diff_rel_df <- rbind(pred_diff_rel_df,data.frame(p_thresh = as.numeric(as.character(comb$thresh[i])),
                                              node_size = as.numeric(comb$node_size[i]),
                                              diff = mean(curDiff),
                                              sd = sd(curDiff),
                                              low_CI = quantile(curDiff,probs = c(0.025)),
                                              high_CI = quantile(curDiff,probs = c(0.975))))
}

pred_diff_abs_df <- data.frame(p_thresh=numeric(),node_size=numeric(),diff=numeric(),sd=numeric())
for(i in 1:length(pred_diff_abs)){
  curDiff <- as.numeric(pred_diff_abs[[i]])
  pred_diff_abs_df <- rbind(pred_diff_abs_df,data.frame(p_thresh = as.numeric(as.character(comb$thresh[i])),
                                                node_size = as.numeric(comb$node_size[i]),
                                                diff = mean(curDiff),
                                                sd = sd(curDiff),
                                                low_CI = quantile(curDiff,probs = c(0.025)),
                                                high_CI = quantile(curDiff,probs = c(0.975))))
}


pred_diff_relpop_df <- data.frame(p_thresh=numeric(),node_size=numeric(),diff=numeric(),sd=numeric())
for(i in 1:length(pred_diff_rel_pop)){
  curDiff <- as.numeric(pred_diff_rel_pop[[i]])
  pred_diff_relpop_df <- rbind(pred_diff_relpop_df,data.frame(p_thresh = as.numeric(as.character(comb$thresh[i])),
                                                        node_size = as.numeric(comb$node_size[i]),
                                                        diff = mean(curDiff),
                                                        sd = sd(curDiff),
                                                        low_CI = quantile(curDiff,probs = c(0.025)),
                                                        high_CI = quantile(curDiff,probs = c(0.975))))
}
library(ggplot2)
pd=position_dodge(0.05)
p <- ggplot(pred_diff_relpop_df,aes(x=-1*log10(p_thresh),y=diff,colour=factor(node_size))) + 
  geom_errorbar(aes(ymin=low_CI,ymax=high_CI),width=.1,position=pd) + 
  geom_line(position=pd)+
  geom_point(position=pd)+
  xlab(expression(paste("Interaction Threshold,","-log"[10],"(p-value)"))) + 
  ylab('Relative Difference in RMSE\n between True and Permuted Phenotype') + 
  labs(colour='Minimum\nNode Size')+
  scale_y_continuous(breaks=seq(-10,0,2),limits=c(-10,0))+
  scale_x_continuous(breaks=seq(3.9,5.4,0.2),limits=c(3.9,5.35)) +
  theme(text = element_text(size=14,family = 'Myriad Pro'))


pd=position_dodge(0)
pred_err_abs_df <- data.frame(p_thresh=as.numeric(as.character(comb$thresh)),node_size=comb$node_size,rmse=sapply(pred_error_abs,function(x) x$nonPermRMSE))
p2 <- ggplot(pred_err_abs_df,aes(x=-1*log10(p_thresh),y=rmse,colour=factor(node_size))) + 
  geom_line(position=pd)+
  geom_point(position=pd)+
  xlab(expression(paste("Interaction Threshold, ","-log"[10],"(p-value)"))) + 
  ylab('Root Mean Standard Error') + 
  labs(colour='Minimum\nNode Size')+
  theme(text = element_text(size=14))

save.image(file='~/bsu_scratch/LDL_Project_Data/Random_Forest/rs12916_rs17238484_rs5909_rs2303152_rs10066707_rs2006760_LDLdirect/pred_err_comparison_plot.RData')

# library(ggplot2)
# p1 <- ggplot(pred_diff_rel_df, aes(x=factor(node_size), y=diff)) +
#   geom_point()+
#   geom_errorbar(aes(ymin=low_CI, ymax=high_CI), width=.2,
#                 position=position_dodge(0.05))+ facet_grid(~factor(p_thresh)) + 
#   labs(y='Mean Relative Difference in Prediction Error',x = 'Minimum Node Size') + theme(text = element_text(size=14))
# 
# p2 <- ggplot(pred_diff_relpop_df, aes(x=factor(node_size), y=diff)) +
#   geom_point()+
#   geom_errorbar(aes(ymin=low_CI, ymax=high_CI), width=.2,
#                 position=position_dodge(0.05))+ facet_grid(~factor(p_thresh))  + 
#   labs(y='Mean Population Relative Difference in Prediction Error',x = 'Minimum Node Size') + theme(text = element_text(size=14)) 
# 
# 
# p3 <- ggplot(pred_diff_abs_df, aes(x=factor(node_size), y=diff)) +
#   geom_point()+
#   geom_errorbar(aes(ymin=low_CI, ymax=high_CI), width=.2,
#                 position=position_dodge(0.05))+ facet_grid(~factor(p_thresh))  + 
#   labs(y='Mean Absolute Difference in Prediction Error',x = 'Min Node Size') + theme(text = element_text(size=14))
