source('~/MRC_BSU_Internship/SNP_Marginal_Effect/pairwise_ld.R')
source('~/MRC_BSU_Internship/Target_Score/gene_phenotype_score.R')
GetLDMatrix <- function(rsids){
  LDMatrix <- matrix(NA,nrow = length(rsids),ncol = length(rsids))
  for(i in 1:nrow(LDMatrix)){
    for(j in 1:i){
      if(i==j){
        LDMatrix[i,j] <- 1
      }else{
        LD <- GetPairwiseLD(rsids[i],rsids[j])
        LDMatrix[i,j] <- ifelse(length(LD)==0,NA,LD)
      }
    }
  }
  LDMatrix[upper.tri(LDMatrix,diag = F)] <- LDMatrix[lower.tri(LDMatrix,diag = F)]
  rownames(LDMatrix) <- rsids
  colnames(LDMatrix) <- rsids
  return(LDMatrix)
}
FindClosestNeighbourLD <- function(ld_matrix,SNP1,SNP2){
  SNP1Neighbours <- ld_matrix[SNP1,]
  SNP2Neighbours <- ld_matrix[SNP2,]
  
  SNP1Neighbours <- names(SNP1Neighbours[SNP1Neighbours > 0.8 & !is.na(SNP1Neighbours)])
  SNP2Neighbours <- names(SNP2Neighbours[SNP2Neighbours > 0.8 & !is.na(SNP2Neighbours)])

  subMatrix <- ld_matrix[SNP1Neighbours,SNP2Neighbours]
  return(ifelse(all(is.na(subMatrix)),NA,max(as.vector(subMatrix),na.rm = T)))
}

IterativePruning <- function(gene_name,upstream_dist,downstream_dist,p_thresh,r2_thresh,n_cores){
  betaCoeff <- CollectBetaCoeff(gene_name,upstream_dist,downstream_dist,0.05,0.5,n_cores)
  betaCoeff <- dplyr::filter(betaCoeff,p < p_thresh)
  if(nrow(betaCoeff) < 2){
    stop('Less than 2 sig SNPs')
  }
  LDMatrix <- GetLDMatrix(betaCoeff$rsid)
  temp <- LDMatrix
  for(i in 1:nrow(LDMatrix)){
    for(j in 1:ncol(LDMatrix)){
      if(is.na(LDMatrix[i,j])){
        LDMatrix[i,j] <- FindClosestNeighbourLD(temp,rownames(LDMatrix)[i],colnames(LDMatrix)[j])
      }
    }
  }
  betaCoeff <- dplyr::arrange(betaCoeff,p)
  inclSet <- betaCoeff$rsid[1]
  for(i in 2:nrow(betaCoeff)){
    ldComparison <- LDMatrix[betaCoeff$rsid[i],inclSet]
    if(any(is.na(ldComparison)) | any(ldComparison > r2_thresh)){
      next
    }
    inclSet <- c(inclSet,betaCoeff$rsid[i])
  }
  return(list(rsid=inclSet,coeff = dplyr::filter(betaCoeff,rsid%in%!!inclSet) %>% {.$coeff}))
}
