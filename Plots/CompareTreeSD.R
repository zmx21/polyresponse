library(pbmcapply)
#Plot comparison of training set and testing set treatment effects (root squared error)
GetSD <- function(path){
  chunks <- c('_1_500','_501_1000','_1001_1500','_1501_2000')
  sdNonPerm <- unlist(lapply(chunks,function(x) readRDS(paste0(path,'sd',x,'.rds'))),recursive = F)
  sdPerm <- unlist(lapply(chunks,function(x) readRDS(paste0(path,'sd_pheno_perm',x,'.rds'))),recursive = F)
  return(list(sdNonPerm=sdNonPerm,sdPerm=sdPerm))
}

resultPath <- '~/bsu_scratch/LDL_Project_Data_Aug2019/Random_Forest_Old/rs12916_rs17238484_rs5909_rs2303152_rs10066707_rs2006760_LDLdirect/'
node_size <- c(10000,20000,30000,40000)
thresh <- c('9e-6','7.5e-6','7.25e-6','7e-6','6.75e-6','6.5e-6','6e-6','5e-6','3e-6')
comb <- expand.grid(node_size,thresh)
thresh <- c('9e-5','7e-5','5e-5','3e-5','1e-5')
node_size <- c(5000,10000,20000,30000,40000)
comb <- rbind(comb,expand.grid(node_size,thresh))


colnames(comb) <- c('node_size','thresh')

# sd_maf_0p05 <- pbmclapply(1:nrow(comb),function(i) GetSD(paste0(resultPath,'0.75_',comb$node_size[i],'_',as.character(comb$thresh[i]),'_5e-2/')),mc.cores = 1)
# saveRDS(sd_maf_0p05,file = '~/bsu_scratch/LDL_Project_Data_Aug2019/Random_Forest_Old/rs12916_rs17238484_rs5909_rs2303152_rs10066707_rs2006760_LDLdirect/tree_sd.rds')
training_sd <- data.frame(p_thresh=numeric(),node_size=numeric(),sd=numeric())
testing_sd <- data.frame(p_thresh=numeric(),node_size=numeric(),sd=numeric(),mean_diff=numeric())
perm_sd <- data.frame(p_thresh=numeric(),node_size=numeric(),sd=numeric(),low_CI=numeric(),high_CI=numeric())

for(i in 1:nrow(comb)){
  curSd <- sd_maf_0p05[[i]]
  curTrainingSd <- sapply(curSd$sdNonPerm,function(x) x$sdTraining)
  curTestingSd <- sapply(curSd$sdNonPerm,function(x) x$sdTesting)
  curPermSd <- curSd$sdPerm

  training_sd <- rbind(training_sd,data.frame(p_thresh = as.numeric(as.character(comb$thresh[i])),
                                              node_size = as.numeric(comb$node_size[i]),
                                              sd = mean(curTrainingSd)))

  testing_sd <- rbind(testing_sd,data.frame(p_thresh = as.numeric(as.character(comb$thresh[i])),
                                              node_size = as.numeric(comb$node_size[i]),
                                              sd = mean(curTestingSd),
                                            mean_diff = mean(sapply(1:length(curTestingSd),function(i) sum(curTestingSd[i] - unlist(curPermSd[[i]]))/length(unlist(curPermSd[[i]])))),
                                            low_CI = quantile(sapply(1:length(curTestingSd),function(i) sum(curTestingSd[i] - unlist(curPermSd[[i]]))/length(unlist(curPermSd[[i]]))),probs = 0.05),
                                            high_CI = quantile(sapply(1:length(curTestingSd),function(i) sum(curTestingSd[i] - unlist(curPermSd[[i]]))/length(unlist(curPermSd[[i]]))),probs = 0.95),
                                            mean_p = mean(sapply(1:length(curTestingSd),function(i) sum(curTestingSd[i] < unlist(curPermSd[[i]]))))))
  curPermSd <- do.call(cbind,curPermSd)
  curPermSd <- sapply(1:nrow(curPermSd),function(x) mean(unlist(curPermSd[x,])))
  perm_sd <- rbind(perm_sd,data.frame(p_thresh = as.numeric(as.character(comb$thresh[i])),
                                            node_size = as.numeric(comb$node_size[i]),
                                            sd = mean(curPermSd),
                                            low_CI = quantile(curPermSd,probs = 0.05),
                                            high_CI = quantile(curPermSd,probs = 0.95)))

}
library(ggplot2)
p1 <- ggplot(training_sd,aes(x=-1*log10(p_thresh),y=sd,colour=factor(node_size))) +
  geom_line()+
  geom_point()+
  xlab(expression(paste("Interaction Threshold, ","-log"[10],"(p-value)"))) +
  ylab('Weighted SD') +
  labs(colour='Minimum\nNode Size') +
  theme(text = element_text(size=14)) +
  scale_y_continuous(breaks=seq(0.0,0.36,0.02),limits=c(0.0,0.36))+
  scale_x_continuous(breaks=seq(3.9,5.4,0.2),limits=c(3.9,5.35)) + ggtitle('Training Set SD')

p2 <- ggplot(testing_sd,aes(x=-1*log10(p_thresh),y=sd,colour=factor(node_size))) +
  geom_line()+
  geom_point()+
  xlab(expression(paste("Interaction Threshold, ","-log"[10],"(p-value)"))) +
  ylab('Weighted SD') +
  labs(colour='Minimum\nNode Size') +
  theme(text = element_text(size=14)) +
  scale_y_continuous(breaks=seq(0.0,0.36,0.02),limits=c(0.0,0.36))+
  scale_x_continuous(breaks=seq(3.9,5.4,0.2),limits=c(3.9,5.35))+ ggtitle('Validation Set SD')

pd=position_dodge(0.1)
p3 <- ggplot(perm_sd,aes(x=-1*log10(p_thresh),y=sd,colour=factor(node_size))) +
  geom_errorbar(aes(ymin=low_CI,ymax=high_CI),width=.1,position=pd) +
  geom_line(position = pd)+
  geom_point(position = pd)+
  xlab(expression(paste("Interaction Threshold, ","-log"[10],"(p-value)"))) +
  ylab('Weighted SD') +
  labs(colour='Minimum\nNode Size') +
  theme(text = element_text(size=14)) +
  scale_y_continuous(breaks=seq(0,0.12,0.02),limits=c(0,0.12))+
  scale_x_continuous(breaks=seq(3.9,5.4,0.2),limits=c(3.9,5.35)) + ggtitle('Permuted Validation Set SD')

p4 <- ggplot(testing_sd,aes(x=-1*log10(p_thresh),y=mean_diff,colour=factor(node_size))) +
  geom_line(position = pd)+
  geom_point(position = pd)+
  xlab(expression(paste("Interaction Threshold, ","-log"[10],"(p-value)"))) +
  ylab('Difference in SD') +
  labs(colour='Minimum\nNode Size') +
  theme(text = element_text(size=14)) + ggtitle('SD Difference between True and Permuted Testing Set') 
p5 <- ggplot(testing_sd,aes(x=-1*log10(p_thresh),y=mean_p/1000,colour=factor(node_size))) +
  geom_line(position = pd)+
  geom_point(position = pd)+
  xlab(expression(paste("Interaction Threshold, ","-log"[10],"(p-value)"))) +
  ylab('P-value, Perm_SD > True_SD') +
  labs(colour='Minimum\nNode Size') +
  theme(text = element_text(size=14)) + ggtitle('SD Difference between True and Permuted Testing Set') 

p6 <- ggplot(dplyr::filter(testing_sd,node_size == 30000),aes(x=-1*log10(p_thresh),y=mean_diff,colour=factor(node_size))) +
  geom_errorbar(aes(ymin=low_CI,ymax=high_CI),width=.1,position=pd)+
  geom_line(position = pd)+
  geom_point(position = pd)+
  xlab(expression(paste("Interaction Threshold, ","-log"[10],"(p-value)"))) +
  ylab('Difference in SD') +
  labs(colour='Minimum\nNode Size') +
  theme(text = element_text(size=14)) + ggtitle('SD Difference between True and Permuted Testing Set') 

p7 <- ggplot(dplyr::filter(testing_sd,node_size == 40000),aes(x=-1*log10(p_thresh),y=mean_diff,colour=factor(node_size))) +
  geom_errorbar(aes(ymin=low_CI,ymax=high_CI),width=.1,position=pd)+
  geom_line(position = pd)+
  geom_point(position = pd)+
  xlab(expression(paste("Interaction Threshold, ","-log"[10],"(p-value)"))) +
  ylab('Difference in SD') +
  labs(colour='Minimum\nNode Size') +
  theme(text = element_text(size=14)) + ggtitle('SD Difference between True and Permuted Testing Set') 

library(ggpubr)
#ggarrange(p2,p3,ncol=2,nrow=1,common.legend = T)
