library(pbmcapply)
#Plot comparison of training set and testing set treatment effects (root squared error)
GetSD <- function(path){
  chunks <- c('_1_500','_501_1000','_1001_1500','_1501_2000')
  sdNonPerm <- unlist(lapply(chunks,function(x) readRDS(paste0(path,'sd',x,'.rds'))),recursive = F)
  sdPerm <- unlist(lapply(chunks,function(x) readRDS(paste0(path,'sd_pheno_perm',x,'.rds'))),recursive = F)
  return(list(sdNonPerm=sdNonPerm,sdPerm=sdPerm))
}

resultPath <- '~/bsu_scratch/LDL_Project_Data/Random_Forest/rs12916_rs17238484_rs5909_rs2303152_rs10066707_rs2006760_LDLdirect/'
node_size <- c(2500,5000,10000,20000,30000,40000)
thresh <- c('5e-6','1e-5','3e-5','5e-5','7e-5','1e-4')
comb <- expand.grid(node_size,thresh)
colnames(comb) <- c('node_size','thresh')

sd <- pbmclapply(1:nrow(comb),function(i) GetSD(paste0(resultPath,'0.75_',comb$node_size[i],'_',as.character(comb$thresh[i]),'/')),mc.cores = 12)
saveRDS(sd,file = '~/bsu_scratch/LDL_Project_Data/Random_Forest/rs12916_rs17238484_rs5909_rs2303152_rs10066707_rs2006760_LDLdirect/tree_sd.rds')
training_sd <- data.frame(p_thresh=numeric(),node_size=numeric(),sd=numeric())
testing_sd <- data.frame(p_thresh=numeric(),node_size=numeric(),sd=numeric())
perm_sd <- data.frame(p_thresh=numeric(),node_size=numeric(),sd=numeric(),low_CI=numeric(),high_CI=numeric())

for(i in 1:nrow(comb)){
  curSd <- sd[[i]]
  curTrainingSd <- sapply(curSd$sdNonPerm,function(x) x$sdTraining)
  curTestingSd <- sapply(curSd$sdNonPerm,function(x) x$sdTesting)
  curPermSd <- curSd$sdPerm
  curPermSd <- do.call(cbind,curPermSd)
  curPermSd <- sapply(1:nrow(curPermSd),function(x) mean(unlist(curPermSd[x,])))
  
  training_sd <- rbind(training_sd,data.frame(p_thresh = as.numeric(as.character(comb$thresh[i])),
                                              node_size = as.numeric(comb$node_size[i]),
                                              sd = mean(curTrainingSd)))
  
  testing_sd <- rbind(testing_sd,data.frame(p_thresh = as.numeric(as.character(comb$thresh[i])),
                                              node_size = as.numeric(comb$node_size[i]),
                                              sd = mean(curTestingSd)))
  
  perm_sd <- rbind(perm_sd,data.frame(p_thresh = as.numeric(as.character(comb$thresh[i])),
                                            node_size = as.numeric(comb$node_size[i]),
                                            sd = mean(curPermSd),
                                            low_CI = quantile(curPermSd,probs = 0.025),
                                            high_CI = quantile(curPermSd,probs = 0.975)))
  
}
library(ggplot2)
p1 <- ggplot(training_sd,aes(x=-1*log10(p_thresh),y=sd,colour=factor(node_size))) +
  geom_line()+
  geom_point()+
  xlab(expression(paste("Interaction Threshold, ","-log"[10],"(p-value)"))) +
  ylab('Weighted SD') +
  labs(colour='Minimum\nNode Size') +
  theme(text = element_text(size=14)) + 
  scale_y_continuous(breaks=seq(0,0.4,0.02),limits=c(0,0.4))+
  scale_x_continuous(breaks=seq(3.9,5.4,0.2),limits=c(3.9,5.35)) + ggtitle('Training Set SD')

p2 <- ggplot(testing_sd,aes(x=-1*log10(p_thresh),y=sd,colour=factor(node_size))) +
  geom_line()+
  geom_point()+
  xlab(expression(paste("Interaction Threshold, ","-log"[10],"(p-value)"))) +
  ylab('Weighted SD') +
  labs(colour='Minimum\nNode Size') +
  theme(text = element_text(size=14)) + 
  scale_y_continuous(breaks=seq(0,0.4,0.02),limits=c(0,0.4))+
  scale_x_continuous(breaks=seq(3.9,5.4,0.2),limits=c(3.9,5.35)) + ggtitle('Validation Set SD')

pd=position_dodge(0.1)
p3 <- ggplot(perm_sd,aes(x=-1*log10(p_thresh),y=sd,colour=factor(node_size))) +
  geom_errorbar(aes(ymin=low_CI,ymax=high_CI),width=.1,position=pd) +
  geom_line(position = pd)+
  geom_point(position = pd)+
  xlab(expression(paste("Interaction Threshold, ","-log"[10],"(p-value)"))) +
  ylab('Weighted SD') +
  labs(colour='Minimum\nNode Size') +
  theme(text = element_text(size=14)) + 
  scale_y_continuous(breaks=seq(0,0.4,0.02),limits=c(0,0.4))+
  scale_x_continuous(breaks=seq(3.9,5.4,0.2),limits=c(3.9,5.35)) + ggtitle('Permuted Validation Set SD')
library(ggpubr)
ggarrange(p1,p2,p3,ncol=3,nrow=1,common.legend = T)
