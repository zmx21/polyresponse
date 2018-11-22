library(pbmcapply)
PredictBetaPerm <- function(path,n_cores,chunks){
  allTrees <- dir(paste0(path,'prediction_betas_perm/'))
  allTrees <- allTrees[sapply(allTrees,function(x) grepl(pattern = 'tree',x=x))]
  
  sampleNodeAssignment <- readRDS(paste0(path,'prediction_betas_perm/',allTrees[1]))$nodeAssignment
  sumBetas <- matrix(0,nrow = nrow(sampleNodeAssignment),ncol = ncol(sampleNodeAssignment))
  pb <- progressBar()
  for(i in chunks){
    setTxtProgressBar(pb,i/length(chunks))
    curTree <- readRDS(paste0(path,'prediction_betas_perm/',allTrees[i]))
    curBeta <- mclapply(1:nrow(sumBetas),function(i){
      curTree$testingMainEffect[i,as.character(curTree$nodeAssignment[i,])]
    },mc.cores = n_cores)
    curBeta <- do.call(rbind,curBeta)
    sumBetas <- sumBetas + curBeta
  }
  close(pb)
  return(sumBetas)
}
resultPath <- '~/bsu_scratch/Random_Forest/rs3821843_rs7340705_rs113210396_rs312487_rs11719824_rs3774530_rs3821856_sbp/'
args <- (commandArgs(T))
node_size <- as.numeric(args[[1]])
p_thresh <- args[[2]]
n_cores <- as.numeric(args[[3]])
path <- paste0(resultPath,'0.75_',node_size,'_',p_thresh,'/')

permResult <- PredictBetaPerm(path,n_cores,1:5000)
saveRDS(permResult,file = paste0(path,'prediction_betas_perm/','mean_beta.rds'))