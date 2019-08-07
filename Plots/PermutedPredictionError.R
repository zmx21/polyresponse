library(pbmcapply)
library(dbplyr)
CalcPredErrorAbs <- function(path,err_type){
  print(path)
  permErrorFiles <- sort(dir(path,pattern = 'beta_err_pheno_perm'))
  nonPermErrorFiles <- sort(setdiff(dir(path,pattern = 'beta_err'),permErrorFiles))
  nonPermError <- unlist(lapply(nonPermErrorFiles,function(x) readRDS(paste0(path,x))),recursive=F)
  
  if(length(permErrorFiles) > 1){
    chunks <- lapply(permErrorFiles,function(x) strsplit(x,'_')[[1]])
    chunks <- lapply(chunks,function(x) sapply(x,function(y) gsub(x = y,pattern = '.rds',replacement = '',fixed = T)))
    chunks <- lapply(chunks,function(x) as.numeric(x)[!is.na(as.numeric(x))])
    chunks <- lapply(chunks,function(x) x[1]:x[2])
  }else{
    chunks <- list(c(1:2000))
    permError <- readRDS(paste0(path,'beta_err_pheno_perm.rds'))
  }
  nonPermErr<- 0
  permErr <- rep(0,1000)
  
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
    if(err_type == 'RMSE'){
      permErr <- permErr + sapply(curPerm,function(x) sum(x^2))
      nonPermErr <-nonPermErr + sum(curNonPerm^2)
    }else if(err_type == 'MAE'){
      permErr <- permErr + sapply(curPerm,function(x) sum(abs(x)))
      nonPermErr <-nonPermErr + sum(abs(curNonPerm))
    }
    numSubGroup <- numSubGroup + length(curNonPerm)
    
  }    
  if(err_type == 'RMSE'){
    return(list(permErr=sqrt(permErr/numSubGroup),nonPermErr=sqrt(nonPermErr/numSubGroup)))
  }else if(err_type == 'MAE'){
    return(list(permErr=permErr/numSubGroup,nonPermErr=nonPermErr/numSubGroup))
  }
}



resultPath <- '~/bsu_scratch/LDL_Project_Data/Random_Forest/rs12916_rs17238484_rs5909_rs2303152_rs10066707_rs2006760_LDLdirect/'
node_size <- c(1000,2500,5000,10000,20000,30000,40000)
thresh <- c('5e-6','1e-5','3e-5','5e-5','7e-5','1e-4')
comb <- expand.grid(node_size,thresh)
colnames(comb) <- c('node_size','thresh')

pred_error_abs_RMSE <- lapply(1:nrow(comb),function(i) CalcPredErrorAbs(paste0(resultPath,'0.75_',comb$node_size[i],'_',as.character(comb$thresh[i]),'/'),'RMSE'))
pred_error_abs_MAE <- lapply(1:nrow(comb),function(i) CalcPredErrorAbs(paste0(resultPath,'0.75_',comb$node_size[i],'_',as.character(comb$thresh[i]),'/'),'MAE'))


pred_error_abs_RMSE_non_perm_df <- data.frame(p_thresh=numeric(),node_size=numeric(),diff=numeric(),sd=numeric())
for(i in 1:length(pred_error_abs_RMSE)){
  curDiff <- as.numeric(pred_error_abs_RMSE[[i]]$nonPermErr)
  pred_error_abs_RMSE_non_perm_df <- rbind(pred_error_abs_RMSE_non_perm_df,data.frame(p_thresh = as.numeric(as.character(comb$thresh[i])),
                                                                    node_size = as.numeric(comb$node_size[i]),
                                                                    diff = mean(curDiff),
                                                                    sd = sd(curDiff),
                                                                    low_CI = quantile(curDiff,probs = c(0.025)),
                                                                    high_CI = quantile(curDiff,probs = c(0.975))))
}
pred_error_abs_RMSE_perm_df <- data.frame(p_thresh=numeric(),node_size=numeric(),diff=numeric(),sd=numeric())
for(i in 1:length(pred_error_abs_RMSE)){
  curDiff <- as.numeric(pred_error_abs_RMSE[[i]]$permErr)
  pred_error_abs_RMSE_perm_df <- rbind(pred_error_abs_RMSE_perm_df,data.frame(p_thresh = as.numeric(as.character(comb$thresh[i])),
                                                                                      node_size = as.numeric(comb$node_size[i]),
                                                                                      diff = mean(curDiff),
                                                                                      sd = sd(curDiff),
                                                                                      low_CI = quantile(curDiff,probs = c(0.025)),
                                                                                      high_CI = quantile(curDiff,probs = c(0.975))))
}

pred_error_abs_RMSE_df <- data.frame(p_thresh=numeric(),node_size=numeric(),diff=numeric(),sd=numeric())
for(i in 1:length(pred_error_abs_RMSE)){
  curDiff <- as.numeric(pred_error_abs_RMSE[[i]]$nonPermErr) - as.numeric(pred_error_abs_RMSE[[i]]$permErr)
  pred_error_abs_RMSE_df <- rbind(pred_error_abs_RMSE_df,data.frame(p_thresh = as.numeric(as.character(comb$thresh[i])),
                                              node_size = as.numeric(comb$node_size[i]),
                                              diff = mean(curDiff),
                                              sd = sd(curDiff),
                                              low_CI = quantile(curDiff,probs = c(0.025)),
                                              high_CI = quantile(curDiff,probs = c(0.975))))
}

pred_error_rel_non_perm_RMSE_df <- data.frame(p_thresh=numeric(),node_size=numeric(),diff=numeric(),sd=numeric())
for(i in 1:length(pred_error_abs_RMSE)){
  curDiff <- (as.numeric(pred_error_abs_RMSE[[i]]$nonPermErr) - as.numeric(pred_error_abs_RMSE[[i]]$permErr)) / as.numeric(pred_error_abs_RMSE[[i]]$nonPermErr)
  pred_error_rel_non_perm_RMSE_df <- rbind(pred_error_rel_non_perm_RMSE_df,data.frame(p_thresh = as.numeric(as.character(comb$thresh[i])),
                                                                    node_size = as.numeric(comb$node_size[i]),
                                                                    diff = mean(curDiff),
                                                                    sd = sd(curDiff),
                                                                    low_CI = quantile(curDiff,probs = c(0.025)),
                                                                    high_CI = quantile(curDiff,probs = c(0.975))))
}

pred_error_rel_perm_RMSE_df <- data.frame(p_thresh=numeric(),node_size=numeric(),diff=numeric(),sd=numeric())
for(i in 1:length(pred_error_abs_RMSE)){
  curDiff <- (as.numeric(pred_error_abs_RMSE[[i]]$nonPermErr) - as.numeric(pred_error_abs_RMSE[[i]]$permErr)) / as.numeric(pred_error_abs_RMSE[[i]]$permErr)
  pred_error_rel_perm_RMSE_df <- rbind(pred_error_rel_perm_RMSE_df,data.frame(p_thresh = as.numeric(as.character(comb$thresh[i])),
                                                                                      node_size = as.numeric(comb$node_size[i]),
                                                                                      diff = mean(curDiff),
                                                                                      sd = sd(curDiff),
                                                                                      low_CI = quantile(curDiff,probs = c(0.025)),
                                                                                      high_CI = quantile(curDiff,probs = c(0.975))))
}


library(ggplot2)
library(ggpubr)
pd=position_dodge(0.1)
p1 <- ggplot(pred_error_abs_RMSE_non_perm_df,aes(x=-1*log10(p_thresh),y=diff,colour=factor(node_size))) +
  geom_errorbar(aes(ymin=low_CI,ymax=high_CI),width=.1,position=pd) +
  geom_line(position=pd)+
  geom_point(position=pd)+
  xlab(expression(paste("Interaction Threshold, ","-log"[10],"(p-value)"))) +
  ylab('RMSE of True Data') +
  labs(colour='Minimum\nNode Size') +
  theme(text = element_text(size=14)) + 
  scale_y_continuous(breaks=seq(0,0.8,0.1),limits=c(0,0.8))+
  scale_x_continuous(breaks=seq(3.9,5.4,0.2),limits=c(3.9,5.35)) 
  
p2 <- ggplot(pred_error_abs_RMSE_perm_df,aes(x=-1*log10(p_thresh),y=diff,colour=factor(node_size))) +
  geom_errorbar(aes(ymin=low_CI,ymax=high_CI),width=.1,position=pd) +
  geom_line(position=pd)+
  geom_point(position=pd)+
  xlab(expression(paste("Interaction Threshold, ","-log"[10],"(p-value)"))) +
  ylab('RMSE of Permuted Data') +
  labs(colour='Minimum\nNode Size') +
  theme(text = element_text(size=14))+
  scale_y_continuous(breaks=seq(0,0.8,0.1),limits=c(0,0.8))+
  scale_x_continuous(breaks=seq(3.9,5.4,0.2),limits=c(3.9,5.35)) 

p3 <- ggplot(pred_error_abs_RMSE_df,aes(x=-1*log10(p_thresh),y=diff,colour=factor(node_size))) +
  geom_errorbar(aes(ymin=low_CI,ymax=high_CI),width=.1,position=pd) +
  geom_line(position=pd)+
  geom_point(position=pd)+
  xlab(expression(paste0("Interaction Threshold, ","-log"[10],"(p-value)"))) +
  ylab('Absolute Difference in RMSE\n between True and Permuted') +
  labs(colour='Minimum\nNode Size') +
  scale_y_continuous(breaks=seq(-0.35,-0.05,0.1),limits=c(-0.35,-0.05))+
  scale_x_continuous(breaks=seq(3.9,5.4,0.2),limits=c(3.9,5.35)) +
  theme_pubr(base_size = 14,base_family = 'Myriad Pro',border = T) #+ scale_color_brewer(type='qual',palette = 'Dark2')

p4 <- ggplot(pred_error_rel_non_perm_RMSE_df,aes(x=-1*log10(p_thresh),y=diff,colour=factor(node_size))) +
  geom_errorbar(aes(ymin=low_CI,ymax=high_CI),width=.1,position=pd) +
  geom_line(position=pd)+
  geom_point(position=pd)+
  xlab(expression(paste("Interaction Threshold, ","-log"[10],"(p-value)"))) +
  ylab('Relative Difference in RMSE\n between True and Permuted based on True') +
  labs(colour='Minimum\nNode Size') +
  theme(text = element_text(size=14))


p5 <- ggplot(pred_error_rel_perm_RMSE_df,aes(x=-1*log10(p_thresh),y=diff,colour=factor(node_size))) +
  geom_errorbar(aes(ymin=low_CI,ymax=high_CI),width=.1,position=pd) +
  geom_line(position=pd)+
  geom_point(position=pd)+
  xlab(expression(paste("Interaction Threshold, ","-log"[10],"(p-value)"))) +
  ylab('Relative Difference in RMSE\n between True and Permuted based on Permuted') +
  labs(colour='Minimum\nNode Size') +
  theme(text = element_text(size=14))

save.image(file='~/bsu_scratch/LDL_Project_Data/Random_Forest/rs12916_rs17238484_rs5909_rs2303152_rs10066707_rs2006760_LDLdirect/pred_err_comparison_plot.RData')

# p2 <- ggplot(pred_error_abs_MAE_df,aes(x=-1*log10(p_thresh),y=diff,colour=factor(node_size))) +
#   geom_errorbar(aes(ymin=low_CI,ymax=high_CI),width=.1,position=pd) +
#   geom_line(position=pd)+
#   geom_point(position=pd)+
#   xlab(expression(paste("Interaction Threshold, ","-log"[10],"(p-value)"))) +
#   ylab('Absolute Difference in MAE\n between True and Permuted') +
#   labs(colour='Minimum\nNode Size') +
#   scale_y_continuous(breaks=seq(-0.35,-0.05,0.1),limits=c(-0.35,-0.05))+
#   scale_x_continuous(breaks=seq(3.9,5.4,0.2),limits=c(3.9,5.35)) +
#   theme(text = element_text(size=14))
# 
# p3 <- ggplot(pred_error_rel_RMSE_df,aes(x=-1*log10(p_thresh),y=diff,colour=factor(node_size))) +
#   geom_errorbar(aes(ymin=low_CI,ymax=high_CI),width=.1,position=pd) +
#   geom_line(position=pd)+
#   geom_point(position=pd)+
#   xlab(expression(paste("Interaction Threshold, ","-log"[10],"(p-value)"))) +
#   ylab('Absolute Difference in RMSE\n between True and Permuted') +
#   labs(colour='Minimum\nNode Size') #+
#   #scale_y_continuous(breaks=seq(-0.35,-0.05,0.1),limits=c(-0.35,-0.05))+
#   #scale_x_continuous(breaks=seq(3.9,5.4,0.2),limits=c(3.9,5.35)) +
#   #theme(text = element_text(size=14))

# print('Calc abs diff')
# pred_diff_abs <- lapply(1:nrow(comb),function(i) CalcPredErrorDiffAbs(paste0(resultPath,'0.75_',comb$node_size[i],'_',as.character(comb$thresh[i]),'/')))
# print('Calc rel diff')
# pred_diff_rel <- lapply(1:nrow(comb),function(i) CalcPredErrorDiffRel(paste0(resultPath,'0.75_',comb$node_size[i],'_',as.character(comb$thresh[i]),'/')))
# print('Calc rel pop diff')
# pred_diff_rel_pop <- lapply(1:nrow(comb),function(i) CalcPredErrorDiffRelPop(paste0(resultPath,'0.75_',comb$node_size[i],'_',as.character(comb$thresh[i]),'/')))
# 
# 
# pred_diff_rel_df <- data.frame(p_thresh=numeric(),node_size=numeric(),diff=numeric(),sd=numeric())
# for(i in 1:length(pred_diff_rel)){
#   curDiff <- as.numeric(pred_diff_rel[[i]])
#   pred_diff_rel_df <- rbind(pred_diff_rel_df,data.frame(p_thresh = as.numeric(as.character(comb$thresh[i])),
#                                               node_size = as.numeric(comb$node_size[i]),
#                                               diff = mean(curDiff),
#                                               sd = sd(curDiff),
#                                               low_CI = quantile(curDiff,probs = c(0.025)),
#                                               high_CI = quantile(curDiff,probs = c(0.975))))
# }
# 
# pred_diff_abs_df <- data.frame(p_thresh=numeric(),node_size=numeric(),diff=numeric(),sd=numeric())
# for(i in 1:length(pred_diff_abs)){
#   curDiff <- as.numeric(pred_diff_abs[[i]])
#   pred_diff_abs_df <- rbind(pred_diff_abs_df,data.frame(p_thresh = as.numeric(as.character(comb$thresh[i])),
#                                                 node_size = as.numeric(comb$node_size[i]),
#                                                 diff = mean(curDiff),
#                                                 sd = sd(curDiff),
#                                                 low_CI = quantile(curDiff,probs = c(0.025)),
#                                                 high_CI = quantile(curDiff,probs = c(0.975))))
# }
# 
# 
# pred_diff_relpop_df <- data.frame(p_thresh=numeric(),node_size=numeric(),diff=numeric(),sd=numeric())
# for(i in 1:length(pred_diff_rel_pop)){
#   curDiff <- as.numeric(pred_diff_rel_pop[[i]])
#   pred_diff_relpop_df <- rbind(pred_diff_relpop_df,data.frame(p_thresh = as.numeric(as.character(comb$thresh[i])),
#                                                         node_size = as.numeric(comb$node_size[i]),
#                                                         diff = mean(curDiff),
#                                                         sd = sd(curDiff),
#                                                         low_CI = quantile(curDiff,probs = c(0.025)),
#                                                         high_CI = quantile(curDiff,probs = c(0.975))))
# }
# library(ggplot2)
# pd=position_dodge(0.05)
# p <- ggplot(pred_diff_relpop_df,aes(x=-1*log10(p_thresh),y=diff,colour=factor(node_size))) + 
#   geom_errorbar(aes(ymin=low_CI,ymax=high_CI),width=.1,position=pd) + 
#   geom_line(position=pd)+
#   geom_point(position=pd)+
#   xlab(expression(paste("Interaction Threshold, ","-log"[10],"(p-value)"))) + 
#   ylab('Relative Difference in RMSE\n between True and Permuted') + 
#   labs(colour='Minimum\nNode Size')+
#   scale_y_continuous(breaks=seq(-10,0,2),limits=c(-10,0))+
#   scale_x_continuous(breaks=seq(3.9,5.4,0.2),limits=c(3.9,5.35)) +
#   theme(text = element_text(size=14))
# 
# 
# pd=position_dodge(0)
# pred_err_abs_df <- data.frame(p_thresh=as.numeric(as.character(comb$thresh)),node_size=comb$node_size,rmse=sapply(pred_error_abs,function(x) x$nonPermRMSE))
# p2 <- ggplot(pred_err_abs_df,aes(x=-1*log10(p_thresh),y=rmse,colour=factor(node_size))) + 
#   geom_line(position=pd)+
#   geom_point(position=pd)+
#   xlab(expression(paste("Interaction Threshold, ","-log"[10],"(p-value)"))) + 
#   ylab('Root Mean Standard Error') + 
#   labs(colour='Minimum\nNode Size')+
#   theme(text = element_text(size=14))
# 

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
