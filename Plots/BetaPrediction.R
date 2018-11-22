library(pbmcapply)
PredictBeta <- function(path,n_cores){
  allTrees <- dir(paste0(path,'prediction_betas/'))
  allTrees <- allTrees[sapply(allTrees,function(x) grepl(pattern = 'tree',x=x))]
  
  sumBetas <- rep(0,length(readRDS(paste0(path,'prediction_betas/',allTrees[1]))$nodeAssignment))
  pb <- progressBar()
  for(i in 1:length(allTrees)){
    setTxtProgressBar(pb,i/length(allTrees))
    curTree <- readRDS(paste0(path,'prediction_betas/',allTrees[i]))
    sumBetas <- sumBetas + curTree$trainingMainEffect[as.character(curTree$nodeAssignment)]
  }
  close(pb)
  return(sumBetas)
}
resultPath <- '~/bsu_scratch/Random_Forest/rs3821843_rs7340705_rs113210396_rs312487_rs11719824_rs3774530_rs3821856_sbp/'
node_size <- 10000
p_thresh <- '3e-5'
n_cores <- 16
path <- paste0(resultPath,'0.75_',node_size,'_',p_thresh,'/')

nonPermResult <- PredictBeta(path,16)