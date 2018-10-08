source('~/MRC_BSU_Internship/Target_Score/gene_phenotype_score.R')
GetLDMatrix <- function(rsids){
  #Add seperator to rsids
  rsids <- paste(rsids,collapse = '%0A')
  #Get result from LDLink tool
  result <- system(command = paste0("curl -k -X GET 'https://analysistools.nci.nih.gov/LDlink/LDlinkRest/ldmatrix?snps=",rsids,"&pop=ALL&r2_d=r2'"),intern = T,ignore.stderr = T)
  #Convert string to df
  result_matrix <- read.table(text = paste(result,collapse = '\n'),header = T)
  #set rownames
  rownames(result_matrix) <- result_matrix$RS_number
  #Remove column of row names
  result_matrix <- result_matrix[,-1]
  #return matrix
  return(as.matrix(result_matrix))
  
  
  # LDMatrix <- matrix(NA,nrow = length(rsids),ncol = length(rsids))
  # for(i in 1:nrow(LDMatrix)){
  #   for(j in 1:i){
  #     if(i==j){
  #       LDMatrix[i,j] <- 1
  #     }else{
  #       LD <- GetPairwiseLD(rsids[i],rsids[j])
  #       LDMatrix[i,j] <- ifelse(length(LD)==0,NA,LD)
  #     }
  #   }
  # }
  # LDMatrix[upper.tri(LDMatrix,diag = F)] <- LDMatrix[lower.tri(LDMatrix,diag = F)]
  # rownames(LDMatrix) <- rsids
  # colnames(LDMatrix) <- rsids
  # return(LDMatrix)
}
# FindClosestNeighbourLD <- function(ld_matrix,SNP1,SNP2){
#   SNP1Neighbours <- ld_matrix[SNP1,]
#   SNP2Neighbours <- ld_matrix[SNP2,]
#   
#   SNP1Neighbours <- names(SNP1Neighbours[SNP1Neighbours > 0.8 & !is.na(SNP1Neighbours)])
#   SNP2Neighbours <- names(SNP2Neighbours[SNP2Neighbours > 0.8 & !is.na(SNP2Neighbours)])
# 
#   subMatrix <- ld_matrix[SNP1Neighbours,SNP2Neighbours]
#   return(ifelse(all(is.na(subMatrix)),NA,max(as.vector(subMatrix),na.rm = T)))
# }

IterativePruning <- function(gene_name,phenotype,upstream_dist,downstream_dist,p_thresh,r2_thresh,n_cores){
  betaCoeff <- CollectBetaCoeff(gene_name,upstream_dist,downstream_dist,phenotype,0.01,0.5,n_cores)
  betaCoeffFilt <- dplyr::filter(betaCoeff,p < p_thresh)
  if(nrow(betaCoeffFilt) < 2){
    print('Less than 2 sig SNPs')
    if(nrow(betaCoeffFilt) == 0){
      betaCoeffFilt <- dplyr::arrange(betaCoeff,p) %>% {.[1,]}
    }
    return(list(rsid=betaCoeffFilt$rsid,coeff = betaCoeffFilt))
    
  }
  print('Finding LD')
  LDMatrix <- GetLDMatrix(betaCoeffFilt$rsid)
  # temp <- LDMatrix
  # for(i in 1:nrow(LDMatrix)){
  #   for(j in 1:ncol(LDMatrix)){
  #     if(is.na(LDMatrix[i,j])){
  #       LDMatrix[i,j] <- FindClosestNeighbourLD(temp,rownames(LDMatrix)[i],colnames(LDMatrix)[j])
  #     }
  #   }
  # }
  betaCoeffFilt <- dplyr::arrange(betaCoeffFilt,p)
  inclSet <- betaCoeffFilt$rsid[1]
  for(i in 2:nrow(betaCoeffFilt)){
    ldComparison = tryCatch({
      LDMatrix[betaCoeffFilt$rsid[i],inclSet]
    }, warning = function(w) {
      ldComparison <- NA
    }, error = function(e) {
      ldComparison <- NA
    })
    if(any(is.na(ldComparison)) | any(ldComparison > r2_thresh)){
      next
    }
    inclSet <- c(inclSet,betaCoeffFilt$rsid[i])
  }
  return(list(rsid=inclSet,coeff = dplyr::filter(betaCoeffFilt,rsid%in%!!inclSet)))
}
