#Plot comparison of root squared error between true testing set and permuted testing set
CalcPredErrorDiff <- function(path){
  print(path)
  nonPermError <- readRDS(paste0(path,'beta_err.rds'))
  permError <- readRDS(paste0(path,'beta_err_perm.rds'))
  percentDiff <- rep(0,1000)

  numSubGroup <- 0
  for(i in 1:length(nonPermError)){
    curNonPerm <- nonPermError[[i]]
    curPerm <- permError[[i]]
    curDiff <- sapply(curPerm,function(x) sum((curNonPerm - x) / abs(curNonPerm)))

    percentDiff <- percentDiff + curDiff
    numSubGroup <- numSubGroup + length(curNonPerm)
  }
  return(percentDiff/numSubGroup)
}

resultPath <- '~/bsu_scratch/Random_Forest/rs3821843_rs7340705_rs113210396_rs312487_rs11719824_rs3774530_rs3821856_sbp/'
node_size <- seq(10000,40000,10000)
thresh <- c('1e-5','2e-5','3e-5','4e-5')
comb <- expand.grid(node_size,thresh)
comb <- comb[-4,]
colnames(comb) <- c('node_size','thresh')


pred_diff <- lapply(1:nrow(comb),function(i) CalcPredErrorDiff(paste0(resultPath,'0.75_',comb$node_size[i],'_',as.character(comb$thresh[i]),'/')))

pred_diff_df <- data.frame(p_thresh=numeric(),node_size=numeric(),diff=numeric(),sd=numeric())
for(i in 1:length(pred_diff)){
  curDiff <- as.numeric(pred_diff[[i]])
  pred_diff_df <- rbind(pred_diff_df,data.frame(p_thresh = as.numeric(as.character(comb$thresh[i])),
                                              node_size = as.numeric(comb$node_size[i]),
                                              diff = mean(curDiff),
                                              sd = sd(curDiff)))
}
library(ggplot2)
p<- ggplot(pred_diff_df, aes(x=factor(node_size), y=diff)) +
  geom_point()+
  geom_errorbar(aes(ymin=diff-sd, ymax=diff+sd), width=.2,
                position=position_dodge(0.05))+ facet_grid(~factor(p_thresh)) + scale_y_continuous(breaks=seq(-14,0,2)) + labs(y='Relative Difference in Prediction Error',x = 'Min Node Size')
