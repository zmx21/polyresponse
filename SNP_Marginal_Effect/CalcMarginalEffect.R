source('/home/zmx21/MRC_BSU_Internship/Load_Bgen/LoadBgen.R')
LinearFit <- function(dosageVector,responseVector,plot){
  #remove missing phenotype samples
  responseVector <- as.numeric(responseVector)
  dosageVector <- as.numeric(dosageVector)
  missingResponse <- is.nan(responseVector)
  responseVector <- responseVector[!missingResponse]
  dosageVector <- dosageVector[!missingResponse]
  data <- data.frame(Y=responseVector,X=dosageVector)
  fit <- lm('Y~X',data = data)
  if(plot){
    boxplot(Y~X,data = data)
  }
  return(data.frame(beta=summary(fit)$coeff[2,1],p=summary(fit)$coeff[2,4]))
}
CalcMarginalEffect <- function(path,samplesTbl,bgen_file_prefix,phenotype,rsid,plot=F){
  dosageVector <- LoadBgen(path,bgen_file_prefix,rsid)
  betaCoeff <- LinearFit(dosageVector,samplesTbl[,phenotype],plot)
  return(cbind(data.frame(rsid=rsid,betaCoeff)))
}
path <-  '/mrc-bsu/scratch/zmx21/UKB_Data/'
sample_file_prefix <- 'ukbb_eur_all_sbp'
bgen_file_prefix <- 'ukb_imp_chr#_HRConly'
phenotype <- 'sbp'
samplesTbl <- LoadSamples(path,sample_file_prefix)

# CalcMarginalEffect(path,samplesTbl,bgen_file_prefix,phenotype,'rs603424')
PKD2L1_marginal <- do.call(rbind,lapply(c('rs603424','rs55825790'),function(x) CalcMarginalEffect(path,samplesTbl,bgen_file_prefix,phenotype,x)))

SCNN1D_marginal <- do.call(rbind,lapply(c('rs1262894','rs1262896','rs1262895'),function(x) CalcMarginalEffect(path,samplesTbl,bgen_file_prefix,phenotype,x)))