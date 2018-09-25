library(dplyr)
library(data.table)
library(rbgen)
library(pbmcapply)
#Obtain reults from PLINK and manually calculated using R
plinkResults <- data.table::fread('/mrc-bsu/scratch/zmx21/UKB_Data/chr10_random_norm/plink.epi.qt')
manualResults <- data.table::fread('/mrc-bsu/scratch/zmx21/UKB_Data/chr10_random_norm/rs603424_sbp/chr10/chr10.txt')

#Sort both results according to rsID
plinkResults <- plinkResults %>% dplyr::select(rsid=SNP2,p_plink=P,coeff_plink=BETA_INT) 
manualResults <- manualResults %>% dplyr::select(rsid=V1,p_manual=V2,coeff_manual=V3,rsq=V4) 

#Load genotype data to look at imputation information
curWd <- getwd()
setwd('/mrc-bsu/scratch/zmx21/UKB_Data/chr10_random_norm/')
chunkSize=50
allRSIds <- unique(plinkResults$rsid)
rsIDChunks <- split(allRSIds,seq(length(allRSIds)-1)%/%chunkSize)

sumOfDeviationList <- pbmclapply(1:length(rsIDChunks),function(i) {
  currentRSIdChunk <- rsIDChunks[[i]]
  genotype_data <-rbgen::bgen.load(filename = 'chr10_random.bgen',rsids = currentRSIdChunk)
  deviation <- pmin(1 - genotype_data$data[,,1],genotype_data$data[,,1]) + pmin(1 - genotype_data$data[,,2],genotype_data$data[,,2]) + pmin(1 - genotype_data$data[,,3],genotype_data$data[,,3])
  sumOfDeviation <- apply(deviation,1,function(x) sum(x)/3/length(x))
  sumOfDeviation <- data.frame(rsid=rownames(genotype_data$data),sum_of_deviation=sumOfDeviation)
  return(sumOfDeviation)
},mc.cores=10)
allResults<- dplyr::inner_join(plinkResults,manualResults,by=c('rsid'='rsid')) %>% dplyr::inner_join(do.call(rbind,sumOfDeviationList),by=c('rsid'='rsid'))
infoScore <- data.table::fread('snp-stats.txt',skip = 9,header = T) %>% dplyr::select(rsid,impute_info)
allResults <- dplyr::inner_join(allResults,infoScore,by=c('rsid'='rsid'))

saveRDS(allResults,'allResults.rds')
setwd(curWd)

noUncertainty <- dplyr::filter(allResults,sum_of_deviation==0)

library(ggplot2)
p1 <- ggplot(data = allResults) + aes(x=coeff_plink,y=coeff_manual) + geom_point() + xlab('PLINK Interaction Coeff') + ylab('R Interaction Coeff') + coord_fixed() + ylim(-0.5,0.5) + xlim(-0.5,0.5) + ggtitle('Interaction coeff comparison')
p2 <- ggplot(data = allResults) + aes(x=-1*log10(p_plink),y=-1*log10(p_manual)) + geom_point() + xlab('PLINK -log10(p-value)') + ylab('R -log10(p-value)') + coord_fixed() + ylim(0,2.5) + xlim(0,2.5) + ggtitle('Interaction p-value comparison')

p3 <- ggplot(data = noUncertainty) + aes(x=-1*log10(p_plink),y=-1*log10(p_manual)) + geom_point() + xlab('PLINK -log10(p-value)') + ylab('R -log10(p-value)') + coord_fixed() + ylim(0,3) + xlim(0,3) + ggtitle('No imputation uncertainty - pvalue comparison')

p4 <- ggplot(data = noUncertainty) + aes(x=coeff_plink,y=coeff_manual) + geom_point() + xlab('PLINK Interaction Coeff') + ylab('R Interaction Coeff') + coord_fixed()  + ggtitle('No imputation uncertainty - interaction coeff comparison')

allResults <- dplyr::arrange(allResults,sum_of_deviation)

p5 <- ggplot(allResults, aes(x=sum_of_deviation, y=abs(-1*log10(p_plink) + 1*log10(p_manual)))) + 
  geom_point() + xlim(0,0.01)  + xlab('Sum of Uncertainty') + ylab('Absolute Difference of log(p-value)') + 
  ggtitle('Uncertainty vs Difference of p-value')

p6 <- ggplot(allResults, aes(x=impute_info, y=abs(-1*log10(p_plink) + 1*log10(p_manual)))) + 
  geom_point() + xlab('Info Score') + ylab('Absolute Difference of log(p-value)') + 
  ggtitle('Info Score vs Difference of p-value') + xlim(0.9,1)
