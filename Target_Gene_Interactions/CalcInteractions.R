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