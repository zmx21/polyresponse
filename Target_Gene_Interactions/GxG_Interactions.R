phenotype = 'sbp'
targetRS <- 'rs603424'

library(dplyr)
source('LoadBgen.R')
#Calculate allele dosage based on genotype probability
alleleProbMatrix <- genotype_data$data
remove(genotype_data)
dosageMatrix <- matrix(0,nrow = nrow(alleleProbMatrix),ncol = ncol(alleleProbMatrix)) + alleleProbMatrix[,,'g=1'] + 2*alleleProbMatrix[,,'g=2']
remove(alleleProbMatrix)

#Construct Sample vs Phenotype Table
samplePhenoTbl <- dplyr::select(samplesTbl,'samples'='ID_1',phenotype)
remove(samplesTbl)
colnames(dosageMatrix) <- samplePhenoTbl$samples

#Function which calculates significance of interaction between two SNPs
CalcInteractions <- function(dosageSubMatrix,phenotypes){
  #Merge into data frame
  data <- data.frame(SNP1=dosageSubMatrix[1,],SNP2=dosageSubMatrix[2,],Pheno=as.numeric(phenotypes))
  #Linear model, with GxG interaction plus individual main effect
  fit <- lm('Pheno~SNP1*SNP2',data = data)
  #Return significance of interaction term
  return(summary(fit)$coefficients[4,4])
}

nonTargetSnps <- 'rs12572586' #setdiff(rownames(dosageMatrix),targetRS)
interactionWithTarget <- mclapply(nonTargetSnps,
                                function(x) CalcInteractions(rbind(dosageMatrix[x,],
                                                                   dosageMatrix[targetRS,]),
                                                             samplePhenoTbl[,phenotype]),mc.cores=20)
