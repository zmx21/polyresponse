library(phenoscanner)
GetPhenoscannerScore <- function(genename,phenotypes){
  res <- phenoscanner(genequery = genename,pvalue = 1,build = 37)
  res <- dplyr::filter(res$results,trait %in% phenotypes)
  return(res)
}