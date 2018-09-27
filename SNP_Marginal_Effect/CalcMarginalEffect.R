source('/home/zmx21/MRC_BSU_Internship/Load_Bgen/LoadBgen.R')
LinearFit <- function(dosageVector,responseVector,plot){
  #remove missing phenotype samples
  responseVector <- as.numeric(responseVector)
  dosageVector <- as.numeric(dosageVector)
  missingResponse <- is.nan(responseVector)
  responseVector <- responseVector[!missingResponse]
  dosageVector <- dosageVector[!missingResponse]
  
  mdlMat <- cbind('Intercept' = rep(1,length(dosageVector)),'Dosage' = dosageVector)
  fit <- RcppEigen::fastLm(y = responseVector,X = mdlMat)
  return(summary(fit)$coeff[2,1])
}
CalcMarginalEffect <- function(path,samplesTbl,bgen_file_prefix,phenotype,rsids,plot=F){
  dosageMatrix <- as.matrix(LoadBgen(path,bgen_file_prefix,rsids))
  betaCoeff <- apply(dosageMatrix,1,function(x) LinearFit(x,samplesTbl[,phenotype],plot))
  return(betaCoeff)
}
# path <-  '/mrc-bsu/scratch/zmx21/UKB_Data/'
# sample_file_prefix <- 'ukbb_eur_all_sbp'
# bgen_file_prefix <- 'ukb_imp_chr#_HRConly'
# phenotype <- 'sbp'
# samplesTbl <- LoadSamples(path,sample_file_prefix)

# CalcMarginalEffect(path,samplesTbl,bgen_file_prefix,phenotype,'rs603424')
