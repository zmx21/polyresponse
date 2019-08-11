####################################################################################
#Helper functional called by GxG_Interactions.R
####################################################################################

#Calculates significance of interaction between two SNPs, including covariates in the model.
CalcInteractions <- function(dosageSubMatrix,dosageTarget,phenotypes,covariates){
  #Create model matrix, with main and interaction effects of two SNPs. 
  mdlMat <- cbind('Intercept' = rep(1,length(dosageSubMatrix)),'snp' = dosageSubMatrix,'target' = dosageTarget,'int' = dosageSubMatrix * dosageTarget,covariates)
  #Fit linear model and calculate stats
  fit <- RcppEigen::fastLmPure(y = phenotypes,X = mdlMat)
  coeff <- fit$coefficients
  se <- fit$se
  t <- coeff/se
  p <- sapply(t,function(x) 2*pt(abs(x),df=fit$df.residual,lower.tail = F))
  #Form matrix, rows are type of coeff, columns are type of variable.
  fit_result <- rbind(coeff,se,t,p)
  #Return coefficient and significance of interaction term
  return(c(as.vector(fit_result[,-1]),fit$s))
}

# CalcInteractions(as.vector(training_set$dosageMatrix[,'rs195781']),training_set$dosageTarget,training_set$phenotypes,as.matrix(training_set$covariates))

# s1 = summary(lm(formula = paste0('LDL ~ HMGCR*rs195781 + ', paste(colnames(training_set$covariates),collapse = '+')),data = cbind(data.frame(LDL = training_set$phenotypes,HMGCR = training_set$dosageTarget,rs195781=training_set$dosageMatrix[,'rs195781']),training_set$covariates)))
# 
# rs195781=LoadBgen(path = '~/bsu_scratch/LDL_Project_Data/Genotype_Data/',bgen_file_prefix = 'ukb_imp_chr#_HRConly',rsIDs = 'rs195781')
# s2 = summary(lm(formula = paste0('LDL ~ HMGCR*rs195781 + ', paste(colnames(training_set$covariates),collapse = '+')),data = cbind(data.frame(LDL = training_set$phenotypes,HMGCR = training_set$dosageTarget,rs195781 = rs195781[,rownames(training_set$dosageMatrix)]),training_set$covariates)))
